#include <tranalisi.hpp>

int main()
{
  input_file_t input("input.txt");
  const size_t _njacks=input.read<size_t>("NJacks");
  set_njacks(_njacks);
  const size_t T=input.read<size_t>("T");
  const size_t ib=input.read<size_t>("ib");
  const size_t norm=input.read<size_t>("norm");
  
  const index_t idLO(vector<pair<string,size_t>>{{"rbw",2},{"rfw",2}});
  const index_t idQED(vector<pair<string,size_t>>{{"rbw0",2},{"rbw1",2},{"rfw0",2},{"rfw1",2}});
  
  const size_t idNoIns1[2]={idLO({0,0}),idLO({0,1})};
  const size_t idNoIns2[2]={idLO({1,1}),idLO({1,0})};
  // const auto& idcc1=idNoIns1;
  // const auto& idcc2=idNoIns2;
  const size_t idDritta1[2]={idQED({0,0,0,0}),idQED({0,0,1,1})};
  const size_t idStorta1[2]={idQED({1,0,1,0}),idQED({1,0,0,1})};
  const size_t idDritta2[2]={idQED({1,1,1,1}),idQED({1,1,0,0})};
  const size_t idStorta2[2]={idQED({0,1,0,1}),idQED({0,1,1,0})};
  
  const string regoTag[2]={"TM","OS"};
  const double zaList[3]={0.703,0.714,0.752},zvList[3]={0.608,0.614,0.657},aList[3]={0.0886,0.0815,0.0619};
  const double za=zaList[ib],zv=zvList[ib],a=aList[ib]/0.197;
  // const double c=za/zv;
  
  const vector<double> x=vector_up_to(T*a/2+1e-3,0.0,a);
  djvec_t connectedSlope;
  
  for(size_t rego=0;rego<2;rego++)
    {
      const auto load=
	[&rego,&T](const string& name,const size_t* a,const size_t* b)
	{
	  const auto load=
	    [&name,&rego,&T](const size_t* i)
	    {
	      return read_djvec("jacks/"+name,T,i[rego]);
	    };
	  
	  return
	    (load(a)+load(b)).symmetrized()/2;
	};
      
      const djvec_t LO_P5P5=load("LO_P5P5",idNoIns1,idNoIns2);
      // const djvec_t LO_cc_P5P5=load("LO_cc_P5P5",idNoIns1,idNoIns2);
      const djvec_t QED_P5P5_dritta=load("QED_P5P5",idDritta1,idDritta2)*zv*zv;
      // const djvec_t QED_cc_P5P5=load("QED_cc_P5P5",idcc1,idcc2);
      const djvec_t QED_P5P5_storta=load("QED_P5P5",idStorta1,idStorta2)*za*za;
      
      const djvec_t rat_dritta=QED_P5P5_dritta/LO_P5P5;
      // const djvec_t rat_cc=QED_cc_P5P5/LO_P5P5;
      const djvec_t rat_storta=QED_P5P5_storta/LO_P5P5;
      
      const djvec_t P5P5_eff=effective_mass(LO_P5P5);
      const djack_t amP5=constant_fit(P5P5_eff,18,47,"plots/LO_P5P5_"+regoTag[rego]+".xmg");
      // const djvec_t cc_P5P5_eff=effective_mass(LO_cc_P5P5);
      const djvec_t CORR_dritta=effective_slope(rat_dritta,P5P5_eff,T/2)/a;
      // const djvec_t CORR_cc=effective_slope(rat_cc,P5P5_eff,T/2);
      const djvec_t CORR_storta=effective_slope(rat_storta,P5P5_eff,T/2)/a;
      if(rego==0)
	connectedSlope=CORR_storta;
      
      const djack_t mP5=amP5/a;
      if(rego==0)
	cout<<"mPi: "<<mP5.ave_err()<<endl;
      
      // cc_P5P5_eff.ave_err().write("plots/LO_cc_P5P5_"+regoTag[rego]+".xmg");
      CORR_dritta.ave_err().write("plots/corr_dritta_"+regoTag[rego]+".xmg");
      // CORR_cc.ave_err().write("plots/corr_cc_"+regoTag[rego]+".xmg");
      CORR_storta.ave_err().write("plots/corr_storta_"+regoTag[rego]+".xmg");
      
      const djvec_t viol=CORR_storta/CORR_dritta-1;
      viol.ave_err().write(x,"plots/viol_"+regoTag[rego]+".xmg");
      
      const djvec_t viol2=forward_derivative(QED_P5P5_storta)/forward_derivative(QED_P5P5_dritta)-1;
      viol2.ave_err().write("plots/viol2_"+regoTag[rego]+".xmg");
      
      // const djvec_t viol_cc=CORR_cc/CORR_dritta/zv/zv-1;
      // viol_cc.ave_err().write("plots/viol_cc_"+regoTag[rego]+".xmg");
      
      // const djvec_t viol2_cc=forward_derivative(QED_cc_P5P5)/forward_derivative(QED_P5P5_dritta)/zv/zv-1;
      // viol2_cc.ave_err().write("plots/viol2_cc_"+regoTag[rego]+".xmg");
      
      CORR_storta.bin_write("CORR_storta_"+regoTag[rego]+".dat");
    }
  
  if(file_exists("jacks/handcuffs"))
    {
      const djvec_t LO_P5P5=read_djvec("jacks/LO_handcuffsRun",T,0).symmetrized()*norm;
      const djvec_t handcuffs=read_djvec("jacks/handcuffs",T).symmetrized()*za*za*norm*norm;
      handcuffs.ave_err().write("plots/handcuffs.xmg");
      const djvec_t handcuffsRat=handcuffs/LO_P5P5;
      handcuffsRat.ave_err().write("plots/handcuffsRat.xmg");
      const djvec_t handcuffsEffSlope=effective_slope(handcuffsRat,effective_mass(LO_P5P5),T/2)/a;
      handcuffsEffSlope.ave_err().write("plots/handcuffsEffSlope.xmg");
      
      const djvec_t totCorr=connectedSlope-handcuffsEffSlope;
      const djack_t tot=constant_fit(totCorr,12,T/2,"plots/totCorr.xmg");
      cout<<"Tot corr: "<<tot.ave_err()<<endl;
    }
  
  return 0;
}
