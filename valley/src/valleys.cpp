
#include <Rcpp.h>
#include <queue>
#include <vector>
#include <functional>
#include <exception>

using namespace Rcpp;



template<class Tx,class Ty> std::pair<Tx,Ty> make_sorted_pair(Tx x,Ty y) {return (x<y)?std::make_pair(x,y):std::make_pair(y,x);}

//' Find valleys of a graph
//' @param x numeric vector of node intensities to consider
//' @param edges integer matrix of 2 columns defining undirected edges of the graph
//' @return an integer vector of size length(x) pointing to the valley of each node
//' @export
//' @author Julien Prados
//' @seealso valleys
//[[Rcpp::export]]
IntegerVector valleysG(NumericVector x,IntegerMatrix edges) {
  // Check parameters
  if (is_true(any(is_na(x)))) throw Rcpp::exception("NA not allowed in x");
  if (edges.ncol()!=2) throw Rcpp::exception("edges must be an integer matrix of 2 columns");
  if (is_true(any(is_na(edges)))) throw Rcpp::exception("NA note allowed in edges");
  if (min(edges)<1 || max(edges)>x.length()) throw Rcpp::exception("edges contains out of range integers");


  // Sort edges in decreasing order considering their lowest value first
  std::vector<int> edges_order(edges.nrow());
  std::iota(edges_order.begin(),edges_order.end(),0);
  std::sort(edges_order.begin(),edges_order.end(),[&x,&edges](int a,int b) {return make_sorted_pair(x[edges(a,0)-1],x[edges(a,1)-1]) > make_sorted_pair(x[edges(b,0)-1],x[edges(b,1)-1]);});


  // Iterate over sorted collection of edges to cluster nodes and compute the valleys
  std::vector<int> top(x.length()); // path to the highest point in the cluster
  std::iota(top.begin(),top.end(),0);
  IntegerVector val(x.length(),NA_INTEGER); // index of the valley of each point
  for(auto i:edges_order) {
    auto e=std::make_pair(edges(i,0)-1,edges(i,1)-1);
    
    // Update cluster membership of the 2 edge's nodes
    while(top[e.first]!=top[top[e.first]]) top[e.first]=top[top[e.first]];
    while(top[e.second]!=top[top[e.second]]) top[e.second]=top[top[e.second]];
    if (top[e.first]==top[e.second]) continue;
    
    // Merge clusters
    if (x[top[e.first]]>x[top[e.second]]) {
      if (IntegerVector::is_na(val[top[e.second]])) val[top[e.second]] = (x[e.first]<x[e.second]?e.first:e.second);
      top[top[e.second]] = top[e.first];
    } else {
      if (IntegerVector::is_na(val[top[e.first]])) val[top[e.first]] = (x[e.first]<x[e.second]?e.first:e.second);
      top[top[e.first]] = top[e.second];
    }
  }

  return val+1;  
}



