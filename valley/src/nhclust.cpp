#include <tuple>
#include <queue>
#include <numeric>
#include <functional>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

struct link_t:std::tuple<double,int,int> {
  link_t(double h,int lhs,int rhs):std::tuple<double,int,int>(h,lhs,rhs) {}
  double& h() {return get<0>(*this);}
  int& lhs() {return get<1>(*this);}
  int& rhs() {return get<2>(*this);}
};


// compute kendall distance between the two given matrix rows
double kendall(NumericMatrix::Row u,NumericMatrix::Row v) {
  // define comparators
  auto cmp1 = [&](size_t i,size_t j){return std::make_pair(u[i],v[i])<std::make_pair(u[j],v[j]);};
  auto cmp2 = [&](size_t i,size_t j){return std::make_pair(v[i],u[i])<std::make_pair(v[j],u[j]);};
  
  // generate a sequence of integer and order it according to cmp1
  std::vector<size_t> o(u.size());
  std::iota(o.begin(),o.end(),0);
  std::sort(o.begin(),o.end(),cmp1);

  // merge sorting of o[] according to cmp2, and count of the number of inverted pairs
  double inv = 0;
  std::vector<size_t> tmp(o.size());
  for (size_t width=1; width<o.size(); width*=2) {
    for(size_t left=0; left<o.size(); left+=2*width) {
      size_t right = std::min(left+width,o.size());
      size_t end = std::min(left+2*width,o.size());
      std::copy(o.begin()+left,o.begin()+end,tmp.begin()+left);
      
      size_t i = left;
      size_t j = right;
      for(size_t k=left;k<end;k++) {
        if (i<right && (j>=end || cmp2(tmp[i],tmp[j]))) {
          o[k] = tmp[i++];
        } else {
          o[k] = tmp[j++];
          inv += right-i;
        }
      }
    }
  }
 
  return inv/(o.size()*(o.size()-1)/2);
}


//' Perform Neighborhood Hierachichal Clustering of matrix rows using centroid linkage and kendall distance
//' @param x numeric matrix to process
//' @param successiveWidth an integer vector splitting rows of x (this vector must sum to nrow(x))
//' @return a data.frame representing the clustering tree: each row represents the merging between 2 clusters, and (nrow(x)-1) merge are necessary for the full clustering.
//'         Field "height" is the kendall distance at which the merging take place, and the data.frame is ordered according to this field.
//'         Fields "(lhs,rhs)" are the indices of merged clusters: a positive integer refers to an already merged cluster, a negative integer refers to a row of x.
//'         Fields "(lb,ub)" are the row range in x spanned by the merged cluster.
//' @seealso \code{\link{nhcentroid}}
//' @export
//' @author Julien Prados
//' @examples
//'   set.seed(999)
//'   x <- matrix(runif(1000),50)
//'   H <- nhclust(x)
//'   with(H,{
//'     plot(range(lb,ub),range(ub-lb+1),type="n")
//'     segments(lb,ub-lb+1,ub,ub-lb+1)
//'   })
// [[Rcpp::export]]
DataFrame nhclust(NumericMatrix x,IntegerVector successiveWidth=IntegerVector(0)) {
  // check parameter "successiveWidth"
  if (successiveWidth.length()==0) successiveWidth.push_back(x.nrow());
  if (accumulate(successiveWidth.begin(),successiveWidth.end(),0) != x.nrow()) throw Rcpp::exception("successiveWidth must sum to nrow(x)");
  if (is_true(any(successiveWidth<0))) throw Rcpp::exception("successiveWidth elements must be >=0");


  // fill priority queue with merging candidates and initialize previous and next cluster informations
  // After this step, q.size() is equal to the expected number of merging step
  std::vector<int> left_element(x.nrow(),-1);
  std::vector<int> right_element(x.nrow(),-1);
  priority_queue<link_t,vector<link_t>,greater<link_t> > q;
  int k=0;
  for(auto n:successiveWidth) {
    for(int i=1;i<n;i++) {
      link_t e(kendall(x.row(k+i-1),x.row(k+i)),k+i-1,k+i);
      left_element[e.rhs()] = e.lhs();
      right_element[e.lhs()] = e.rhs();
      q.push(e);
    }
    k+=n;
  }
  

  // declare an union-find structure as an array of parents, and implement corresponding methods
  std::vector<int> parents(x.nrow()+q.size());
  std::iota(parents.begin(),parents.end(),0);
  auto is_root = [&parents](int gid) {return gid==parents[gid];};
  auto find_root = [&parents](int gid) {
    while (gid != parents[gid]) gid = parents[gid];
    return gid;
  };
  
  
  // a method to retreive the centroid_row of a cluster
  NumericMatrix mu(q.size(),x.ncol());  
  auto centroid_row = [&x,&mu](int gid) {return (gid<x.nrow())?x.row(gid):mu.row(gid-x.nrow());};


  // structure to store boundaries of the clusters and methods to retreive the values as well as cluster sizes
  IntegerVector lb(q.size());
  IntegerVector ub(q.size());
  auto get_lb = [&x,&ub,&lb](int gid) {return (gid<x.nrow())?gid:lb[gid-x.nrow()];};
  auto get_ub = [&x,&ub,&lb](int gid) {return (gid<x.nrow())?gid:ub[gid-x.nrow()];};
  auto cluster_size = [&x,&ub,&lb](int gid) {return (gid<x.nrow())?1:1 + ub[gid-x.nrow()] - lb[gid-x.nrow()];};


  // additional structures to store clusters informations
  NumericVector height(q.size());
  IntegerVector lhs(q.size());
  IntegerVector rhs(q.size());  
  
  
  // iterate throught merging candidates, from most correlated to less correlated
  for(int step = 0;!q.empty();) {
    link_t top = q.top();
    q.pop();
    if (is_root(top.lhs()) && is_root(top.rhs())) {
      // compute boundary information for the new merged cluster
      height[step] = top.h();
      lhs[step] = top.lhs();
      rhs[step] = top.rhs();
      lb[step] = get_lb(top.lhs());
      ub[step] = get_ub(top.rhs());
          
      // compute weighted average of left and right clusters as the centroid of the new merged cluster
      int lsize = cluster_size(top.lhs());
      int rsize = cluster_size(top.rhs());      
      NumericMatrix::Row left_row = centroid_row(top.lhs());
      NumericMatrix::Row right_row = centroid_row(top.rhs());      
      transform(left_row.begin(),left_row.end(),right_row.begin(),mu.row(step).begin(),[=](double x,double y){return (x*lsize + y*rsize)/(lsize+rsize);});
      
      // set the new merged cluster we just create as the parent of left and right clusters
      parents[lhs[step]] = parents[rhs[step]] = parents[lb[step]] = parents[ub[step]] = x.nrow()+step;

      // push new candidate merging the new cluster with previous one
      if (left_element[lb[step]]>=0) {
        top.lhs() = find_root(left_element[lb[step]]);
        top.rhs() = x.nrow() + step;
        top.h() = kendall(centroid_row(top.lhs()),centroid_row(top.rhs()));
        q.push(top);
      }

      // push new candidate merging the new cluster with next one
      if (right_element[ub[step]]>=0) {
        top.lhs() = x.nrow() + step;
        top.rhs() = find_root(right_element[ub[step]]);
        top.h() = kendall(centroid_row(top.lhs()),centroid_row(top.rhs()));
        q.push(top);        
      }
      step++;
    }
  }
  return DataFrame::create(
    Named("height")=height,
    Named("lhs")=ifelse(lhs<x.nrow(),lhs+1,-(lhs-x.nrow()+1)),
    Named("rhs")=ifelse(rhs<x.nrow(),rhs+1,-(rhs-x.nrow()+1)),
    Named("lb")=lb+1,Named("ub")=ub+1
  );
}



//' Compute centroid of an nhclust tree on given sample matrix
//' @param x numeric matrix of the samples to cluster (one sample per column). Number of row in the matrix must be compatible with the provided clustering tree.
//' @param tree a data.frame representing a nhclust tree. This object my be generated with function nhclust(), and eventualy reduced to head rows only.
//' @return a numeric matrix of size nrow(tree) x ncol(x). Each row being the centroid of cluster i, for all samples.
//' @seealso \code{\link{head}}, \code{\link{nhclust}}
//' @export
//' @author Julien Prados
//' @examples
//'   set.seed(999)
//'   x <- matrix(runif(1000),50)
//'   x2 <- matrix(runif(500),50)
//'   H <- nhclust(x)
//'   mu <- nhcentroid(H,x)
//'   mu2 <- nhcentroid(H,x2)
// [[Rcpp::export]]
NumericMatrix nhcentroid(DataFrame tree,NumericMatrix x) {
    IntegerVector lhs = tree["lhs"];
    IntegerVector rhs = tree["rhs"];
    IntegerVector lb = tree["lb"];
    IntegerVector ub = tree["ub"];    
    if (is_true(any(lhs>x.nrow())) || is_true(any(-lhs>tree.nrows()))) throw Rcpp::exception("index out of bound in lhs");
    if (is_true(any(rhs>x.nrow())) || is_true(any(-rhs>tree.nrows()))) throw Rcpp::exception("index out of bound in rhs");
    if (is_true(any(lb>x.nrow())) || is_true(any(lb<0))) throw Rcpp::exception("index out of bound in lb");
    if (is_true(any(ub>x.nrow())) || is_true(any(ub<0))) throw Rcpp::exception("index out of bound in ub");

    NumericMatrix mu(tree.nrows(),x.ncol());
    for(int i=0;i<mu.nrow();i++) {
      NumericMatrix::Row u = (lhs[i]>0)?x.row(lhs[i]-1):mu.row(-lhs[i]-1);
      NumericMatrix::Row v = (rhs[i]>0)?x.row(rhs[i]-1):mu.row(-rhs[i]-1);
      int lsize = (lhs[i]>0)?1:1+ub[-lhs[i]-1]-lb[-lhs[i]-1];
      int rsize = (rhs[i]>0)?1:1+ub[-rhs[i]-1]-lb[-rhs[i]-1];
      transform(u.begin(),u.end(),v.begin(),mu.row(i).begin(),[=](double x,double y){return (x*lsize + y*rsize)/(lsize+rsize);});
    }
    return mu;
}




