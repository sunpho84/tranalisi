#include <tranalisi.hpp>

const int T=48,L=24;
int TH;

int main(){
  set_njacks(15);
  vector<int> nconfs;
  nconfs.push_back(15);
  nconfs.push_back(30);
  nconfs.push_back(75);
  nconfs.push_back(150);
  
  index_t ind({{"ind1",1},{"ind2",1}});
  vector<raw_file_t> fouthiterr(4);
  for(int a=0;a<nconfs.size();a++)
  fouthiterr[a].open(combine("prova%d",a),"w");
  return 0;
}
