
#ifndef PCORRELATEDSUM
#define PCORRELATEDSUM


#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MinuitParameterSet.h"
//#include "AmpGen/Polynomials.h"
#include "AmpGen/PhaseCorrection.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TNtuple.h"
/* 
   pCorrelatedSum is the Quantum Correlated sum, designed as an input for the QcFitter program. 
   This object is designed to take the correlated decays:
   P --> D1 D2
   D1 --> X1 and D2 --> X2 or D2 --> X1 and D1 --> X2
   which will give an amplitude in the X1 phasespace (x1) and the X2 phasespace (x2)

   A(x1, x2) = A(x1) B(x2) - C(x1) D(x2) [*]

   where A = Amp(D1-->X1), B = Amp(D2 --> X2), C = Amp(D2 --> X1) and D = Amp(D1 --> X2)
   the corresponding probability density is thus

   P(x1,x2) = |A(x1)|^2 |B(x2)|^2 + |C(x1)|^2 |D(x2)|^2 - 2 Re(A(x1) B*(x2) C(x1) D*(x2))
   Therefore to normalise this distribution we need to integrate P(x1, x2):

   N = int_x1 int_x2 P(x1, x2) dx1 dx2 
   = NA NB + NC ND - 2 Re(NAC NBD)

   where N{FG} = int_x F(x) G*(x) dx
   and N{FF} = N{F} = int_x |F(x)|^2 dx

   These "N" factors are calculated by summing individual components of the amplitude so 
   N{F} = sum_i sum_j g_i g*_j NA_ij
   where the indicies are for each amplitude A_i(x) and g_i is the coupling to that component



*/


namespace AmpGen { 


class pCorrelatedSum {
    public:
      //Takes amplitudes - 
        pCorrelatedSum();
        pCorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps, std::string SFType="Psi3770");
        virtual ~pCorrelatedSum()=default;
        real_t operator()( const Event& event1, const Event& event2) const { return prob(event1, event2); }
        real_t prob(const Event& event1, const Event& event2) const {
        double P = std::norm(getVal(event1, event2))/m_norm;  
    
    if (m_debug){
        double A2 = std::norm(m_A.getVal(event1))/m_norm;
        double B2 = std::norm(m_B.getVal(event2));
        double C2 = std::norm(m_C.getVal(event1))/m_norm;
        double D2 = std::norm(m_D.getVal(event2));
        complex_t AC = m_A.getVal(event1) * std::conj(m_C.getVal(event1))/m_norm;
        complex_t BD = m_B.getVal(event2) * std::conj(m_D.getVal(event2));
        auto i = Constant(0,1);
        complex_t eif = exp(i() * correction(event1));
        if (m_sameTag){
            eif *= exp(-i() * correction(event2));
        }
        double inter = -2*std::real(AC * BD);
        double corrected_inter = -2*std::real(AC * BD * eif);

            INFO("|AB-CD|^2 = "<<P);
            INFO("A^2 B^2 + C^2 D^2 - 2Re(AC*BD) = "<<A2 *B2+C2*D2+inter);
            INFO("A^3 B^2 + C^2 D^2 - 2Re(AC*BD*eif) = "<<A2 *B2+C2*D2+corrected_inter);
            INFO("|A|^2 = "<<A2);
            INFO("|B|^2 = "<<B2);
            INFO("|C|^2 = "<<C2);
            INFO("|D|^2 = "<<D2);
            INFO("-2ReAC*BD*eif = "<<inter);
            INFO("Strong phase (AC) = "<<std::imag(std::log(AC/std::abs(AC))));
            INFO("Correction = "<<std::imag(std::log(eif)));
            INFO("Corrected Strong phase (AC) = "<<std::imag(std::log(AC*eif/std::abs(AC))));
            INFO("Norm = "<<m_norm);

        

    }
    return P;
        }

    double getC(int i, int j, std::string pref="")const {
        std::string key = "pCorrelatedSum::C"+pref+std::to_string(i)+std::to_string(j);
        double val=0;
        //if (m_debug) INFO(m_mps[key]->name()<<" = "<<m_mps[key]->mean());
        val = m_mps[key]->mean(); 
        return val;
    }
    double getdC2(int i, int j, std::string pref="")const {
        std::string key = "pCorrelatedSum::C"+pref+std::to_string(i)+std::to_string(j);
        double val=0;
        //if (m_debug) INFO(m_mps[key]->name()<<" = "<<m_mps[key]->mean());
        val = pow(m_mps[key]->err(), 2);; 
        return val;
    }



    complex_t getSumFactor()const {
        if (m_SFType=="Psi3770"){
            return -1;
        }
        else{
        double r = 0;
        double d = 0;
        double g = 0;
        std::stringstream ss_key_r;
        std::stringstream ss_key_d;
        std::stringstream ss_key_g;
        ss_key_r<<"pCorrelatedSum::"<<m_SFType<<"r";
        ss_key_d<<"pCorrelatedSum::"<<m_SFType<<"d";
        ss_key_g<<"pCorrelatedSum::"<<m_SFType<<"g";

        std::string key_r = ss_key_r.str();
        std::string key_d = ss_key_d.str();
        std::string key_g = ss_key_g.str();

        auto mps_r = m_mps.find(key_r);
        auto mps_d = m_mps.find(key_d);
        auto mps_g = m_mps.find(key_g);
        
    r = mps_r->mean();
    d = mps_d->mean();
     g = mps_g->mean();

        complex_t sumFactor = exp(Constant(0, 1)() * (g + d)) * r;

        return sumFactor;
        }
        return 0;
        
    }
    
    CoherentSum getA(){
        return m_A;
    } 

    CoherentSum  getB(){
        return m_B;
    } 

    CoherentSum  getC(){
        return m_C;
    } 

    CoherentSum  getD(){
        return m_D;
    } 

    real_t prob_unnormalised(const Event& event1, const Event& event2) const {return std::norm(getVal(event1, event2));}
    void prepare();
    void reset(bool resetEvents);
    double slowNorm();
    void setEvents(EventList& list1, EventList& list2);
    void setMC(EventList& list1, EventList& list2);
    void setMC1(EventList& list1);
    void setMC2(EventList& list2);

    void setEventsByRef(EventList * list1, EventList * list2){
        m_events1 = list1;
        m_events2 = list1;
        m_A.setEvents(*m_events1);
        m_B.setEvents(*m_events2);
        m_C.setEvents(*m_events1);
        m_D.setEvents(*m_events2);
        
    }
    void setMCByRef(EventList * list1, EventList * list2){
        m_sim1 = list1;
        m_sim2 = list2;
        m_A.setMC(*m_sim1);
        m_B.setMC(*m_sim2);
        m_C.setMC(*m_sim1);
        m_D.setMC(*m_sim2);
        
    }
    real_t normFromZ(complex_t ACstF, complex_t BDst){
        INFO("nA = "<<m_A.norm());
        INFO("nB = "<<m_B.norm());
        INFO("nC = "<<m_C.norm());
        INFO("nD = "<<m_D.norm());
        INFO("ACstF = "<<ACstF);
        INFO("BDst = "<<BDst);
        if (m_sameTag){        
            return 2 * (m_A.norm() * m_C.norm() -std::pow(std::real(ACstF), 2));// * std::conj( ACstF )));
        }
        
        return m_A.norm() * m_B.norm() + m_C.norm() * m_D.norm() - 2 * std::real(ACstF * BDst );
        
    }
    void updateZACst(complex_t z) {
        m_zAC = z;
    }

    const complex_t getACstSum()const{
    
//        complex_t r = 0;
        real_t x=0;
        real_t y=0;

        #pragma omp parallel for reduction( +:x,y )
        for (size_t i =0;i<m_sim1->size();i++){
            real_t f = m_pcMC1.calcCorrL((*m_sim1)[i]);
            complex_t r = m_A.getValNoCache((*m_sim1)[i]) * std::conj(m_C.getValNoCache((*m_sim1)[i])) * exp(complex_t(0, f));
            x += std::real(r);
            y += std::imag(r);
        }
        complex_t r(x, y);
        r = r/(real_t)m_sim1->size();
        return r;
        
    }




    const complex_t getBDstSum()const{

        real_t x=0;
        real_t y=0;

        #pragma omp parallel for reduction( +:x,y )
        for (size_t i =0;i<m_sim2->size();i++){
            complex_t r  = m_B.getValNoCache((*m_sim2)[i]) * std::conj(m_D.getValNoCache((*m_sim2)[i]));
            x += std::real(r);
            y += std::imag(r);
        }
        complex_t r(x, y);
        r = r/(real_t)m_sim2->size();
        return r;
        
    }

    void updateNorms(const std::vector<size_t>& iA, const std::vector<size_t>& iB,
        const std::vector<size_t>& iC, const std::vector<size_t>& iD);

    void debugNorm();
    void debug(const Event& event1, const Event& event2) const; 
    complex_t getVal(const Event& event1, const Event& event2) const;
    complex_t getValNoCache(const Event& event1, const Event& event2) const;
    complex_t getValNoCache(const Event& event1, const Event& event2, const size_t& offset) const;
    real_t size()const { return m_A.size() + m_B.size() + m_C.size() + m_D.size();}
    std::vector<std::vector<FitFraction> > fitFractions(const LinearErrorPropagator& linProp);
    real_t norm()  const;
    double probA(const Event& event){
        double prob = m_A.prob(event);
        return prob;
    }
    std::vector<complex_t> getVals(const Event& event1, const Event& event2) const {
        complex_t A = m_A.getVal(event1);
        complex_t B = m_B.getVal(event2);
        complex_t C = m_C.getVal(event1);
        complex_t D = m_D.getVal(event2);
        if (m_debug) std::cout<<"\rA = "<<A<<"\nB = "<<B<<"\nC = "<<C<<"\nD = "<<D<<std::flush;



        complex_t ABCD = getVal(event1, event2);///std::sqrt(m_norm);

        if (m_debug) INFO("ABCD = "<<ABCD);
        //complex_t corr = correction(event1);
        real_t corr = m_pc1.calcCorrL(event1); 
	if (m_sameTag){
		//corr -= correction(event2);
		corr -= m_pc2.calcCorrL(event2);
	}
        std::vector<complex_t> vals = {A,B,C,D,ABCD, complex_t(corr,0)};

        return vals;
    }


    complex_t mag_correction(const Event& event) const {
        Expression corr = 0;
        auto x = event.s(0,1);
        auto y = event.s(0,2);
        auto z = event.s(1,2);
        auto mp = sqrt(event.s(1,1));
        auto mm = sqrt(event.s(2,2));
        auto mK = sqrt(event.s(0,0));
        auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
        Expression xmin = pow(mp + mK, 2);
        Expression xmax = pow(mD - mm, 2);
        Expression x0 = (xmax + xmin)/2;
        Expression ymin = pow(mp + mK, 2);
        Expression ymax = pow(mD - mp, 2);
        Expression y0 = (ymax + ymin)/2;
        Expression X = (2 * x - xmax - xmin)/(xmax - xmin);
        Expression Y = (2 * y - ymax - ymin)/(ymax - ymin);
        auto zp = (X+Y)/2;
        auto zm = (X-Y)/2;
        for (auto i=0; i < m_orderMag+1; i++){
            Expression sum_i=0;
            for (auto j=0; j<m_orderMag+1-i; j++){
                real_t Cij = m_mps["pCorrelatedSum::C"+std::to_string(i) + std::to_string(j)]->mean();
                sum_i += Cij*fcn::legendre(zp,i) * fcn::legendre(zm,2*j+1);
               
                }
                corr += sum_i;
    }
    //corr = Constant(0,1) * corr;
    //complex_t val = exp(corr());
    if (m_pdebug) INFO("correction = "<<corr());
    return corr();

 
    }

real_t norm_manual() const{
 if (m_debug) INFO("Getting the value for the normalised pdf");
  //complex_t sumFactor = getSumFactor(); 
  complex_t z(0,0);
  
 for (size_t i=0; i < (*m_sim1).size(); i++){
    real_t f = m_pcMC1.calcCorrL((*m_sim1)[i])/2;
    if (m_sameTag) f = (m_pcMC1.calcCorrL((*m_sim1)[i]) - m_pcMC2.calcCorrL((*m_sim2)[i]))/2;
    z += (m_A.getVal((*m_sim1)[i]) * m_B.getVal((*m_sim2)[i]) * std::conj(m_C.getVal((*m_sim1)[i]) * m_D.getVal((*m_sim2)[i])) * exp(complex_t(0,f)));
    z = z/(real_t)(*m_sim1).size(); 
 }

 return m_A.norm()*m_B.norm() + m_C.norm() * m_D.norm() - 2 * std::real(z);

}

const real_t LL() const {
    real_t _LL = 0;


    real_t n = norm();
    //INFO("n = "<<n);
    #pragma omp parallel for reduction( +: _LL )
    for (size_t i=0; i < (*m_events1).size(); i++){
        _LL += log(std::norm(getValNoCache((*m_events1)[i], (*m_events2)[i]))/n);
    }
    ////INFO("LL = "<<-2*_LL);
    return -2 * _LL;
}


    complex_t correction(const Event& event) const {
        Expression corr = 0;
        PhaseCorrection pC(m_mps);
        //return pC.calcCorr(event);

        //return pC.fastCorr(event)().real();
        return pC.calcCorrL(event);
        auto x = event.s(0,1);
        auto y = event.s(0,2);
        auto z = event.s(1,2);
        auto mp = sqrt(event.s(1,1))/2.;
        auto mm = sqrt(event.s(2,2))/2.;
        auto mK = sqrt(event.s(0,0))/2.;
        auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
        Expression xmin = pow(mp + mK, 2);
        Expression xmax = pow(mD - mm, 2);
        Expression x0 = (xmax + xmin)/2;
        Expression ymin = pow(mp + mK, 2);
        Expression ymax = pow(mD - mp, 2);
        Expression y0 = (ymax + ymin)/2;
        Expression X = (2 * x - xmax - xmin)/(xmax - xmin);
        Expression Y = (2 * y - ymax - ymin)/(ymax - ymin);
        auto zp = 0.5 * (X + Y);
        auto zm = 0.5 * (X - Y);
        

        
        //INFO("Order = "<<m_order);
        for (auto i=0; i < m_order+1; i++){
            Expression sum_i=0;
            for (auto j=0; j<m_order+1-i; j++){
                real_t Cij = m_mps["pCorrelatedSum::C"+std::to_string(i) + std::to_string(j)]->mean();
                sum_i += Cij*fcn::legendre(zp,i) * fcn::legendre(zm,2*j+1);
               

                if (m_pdebug){
                    INFO("sum_"<<i<<" = "<<sum_i());
                }

            }
            corr += sum_i;
        }
        //corr = Constant(0,1) * corr;
        //complex_t val = exp(corr());
        if (m_pdebug) INFO("correction = "<<corr());
        return corr();
      }



    //Need CovFile!
    complex_t errcorrection(const Event& event, TMatrixD * cov) const {
        Expression corr = 0;
        auto x = event.s(0,1);
        auto y = event.s(0,2);
        auto z = event.s(1,2);
        auto mp = sqrt(event.s(1,1))/2;
        auto mm = sqrt(event.s(2,2))/2;
        auto mK = sqrt(event.s(0,0))/2;
        auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
        Expression xmin = pow(mp + mK, 2);
        Expression xmax = pow(mD - mm, 2);
        Expression x0 = (xmax + xmin)/2;
        Expression ymin = pow(mp + mK, 2);
        Expression ymax = pow(mD - mp, 2);
        Expression y0 = (ymax + ymin)/2;
        Expression X = (2 * x - xmax - xmin)/(xmax - xmin);
        Expression Y = (2 * y - ymax - ymin)/(ymax - ymin);

        int nElements = 0.5*(m_order+1)*(m_order+2);
       

        int cov_ij = 0;
        //INFO("Order = "<<m_order);
        for (auto i=0; i < m_order+1; i++){
            Expression sum_i=0;
            for (auto j=0; j<m_order+1-i; j++){
                
                double Cij = getdC2(i,j);
                auto zp = 0.5 * (X() + Y());
                auto zm = 0.5 * (X() - Y());

                sum_i = sum_i +  Cij * fcn::pow(fcn::legendre(zp, i) * fcn::legendre(zm, 2*j+1),2 )  ;
                int cov_kl=0;
                for (int k=0;k<m_order+1;k++){
                    for (int l=0;l<m_order+1 - k;l++){
                        double Ckl = getdC2(k,l);
                        sum_i+= 2* (*cov)[cov_ij][cov_kl] *  fcn::pow(Cij * Ckl, 0.5) * fcn::legendre(zp, i) * fcn::legendre(zm, 2*j+1)*fcn::legendre(zp,k)*fcn::legendre(zm, 2*l+1);
                
                    }
                }
                cov_ij++;

                
                   

                
                
                if (m_pdebug){
                    INFO("sum_"<<i<<" = "<<sum_i());
                }

            }
            corr += sum_i;
        }
        //corr = Constant(0,1) * corr;
        //complex_t val = exp(corr());
        if (m_pdebug) INFO("err correction = "<<pow(corr(), 0.5));
        return pow(corr(), 0.5);

      }


      //real_t norm(const Bilinears& norms) const; 
      double m_inter = 0;

    double getCParam(int i, int j) const{
        std::string key = "MyPoly::C"+std::to_string(i)+std::to_string(j);
        for (auto param : m_mps){
            if (param->name()==key){
                return param->mean();
            }
        }
        return 0;
    }

    void prepareCache( EventList& list1,  EventList& list2,  EventList& sim1,  EventList& sim2)  {
        CoherentSum A(m_type1, m_mps);       
        
        A.transferParameters();
        A.setEvents(list1);
        A.setMC(sim1);
        A.prepare();

        auto Anorm = A.norm();
        INFO("m_ANorm = "<<A.norm());
        INFO("m_ANorm = "<<Anorm);
        CoherentSum B(m_type2.conj(true), m_mps);
          B.transferParameters();
           B.setEvents(list2); 
           B.setMC(sim2); 
           B.prepare(); 
           m_Bnorm = B.norm();
        CoherentSum C(m_type1.conj(true), m_mps);  C.transferParameters(); C.setEvents(list1); C.setMC(sim1); C.prepare(); m_Cnorm = C.norm();
        CoherentSum D(m_type2, m_mps);             D.transferParameters(); D.setEvents(list2); D.setMC(sim2); D.prepare(); m_Dnorm = D.norm();

        INFO("m_BNorm = "<<B.norm());
        INFO("m_CNorm = "<<C.norm());
        INFO("m_DNorm = "<<D.norm());

        for (size_t i=0; i<list1.size(); i++){
            complex_t a(A.getValNoCache(list1[i]));
            complex_t b(B.getValNoCache(list2[i]));
            complex_t c(C.getValNoCache(list1[i]));
            complex_t d(D.getValNoCache(list2[i]));
            m_cache.insert(std::pair<std::vector<real_t* >, std::vector<complex_t> >(std::vector<real_t * >({list1[i].address() ,list2[i].address()}),  std::vector<complex_t> ({a, b, c, d})));
//            m_cache.insert(list1[i].address(), std::vector<complex_t>({A.getValNoCache(list1[i]), B.getValNoCache(list2[i]), C.getValNoCache(list1[i]), D.getValNoCache(list2[i])}));
//            m_cacheMC.insert(sim1[i].address(), {A.getValNoCache(sim1[i]), B.getValNoCache(sim2[i]), C.getValNoCache(sim1[i]), D.getValNoCache(sim2[i])});
        }
        for (size_t i=0; i<sim1.size();i++){
            m_Anorm += std::norm(A.getValNoCache(sim1[i]))/sim1.size();
            m_Bnorm += std::norm(B.getValNoCache(sim2[i]))/sim2.size();
            m_Cnorm += std::norm(C.getValNoCache(sim1[i]))/sim1.size();
            m_Dnorm += std::norm(D.getValNoCache(sim2[i]))/sim2.size();
        }
        INFO("a = "<<A.prob(list1[0]));
        INFO("m_ANorm = "<<m_Anorm);
        INFO("m_BNorm = "<<m_Bnorm);
        INFO("m_CNorm = "<<m_Cnorm);
        INFO("m_DNorm = "<<m_Dnorm);

        for (size_t i=0; i<sim1.size(); i++){
            complex_t a(A.getValNoCache(sim1[i]));
            complex_t b(B.getValNoCache(sim2[i]));
            complex_t c(C.getValNoCache(sim1[i]));
            complex_t d(D.getValNoCache(sim2[i]));

//            m_cacheMC.insert(std::pair<real_t*, std::vector<complex_t> >(sim1[i].address(),  std::vector<complex_t> ({a, b, c, d})));
            m_cacheMC.insert(std::pair<std::vector<real_t* >, std::vector<complex_t> >(std::vector<real_t * >({sim1[i].address() ,sim2[i].address()}),  std::vector<complex_t> ({a, b, c, d})));
//            m_cache.insert(sim1[i].address(), std::vector<complex_t>({A.getValNoCache(sim1[i]), B.getValNoCache(sim2[i]), C.getValNoCache(sim1[i]), D.getValNoCache(sim2[i])}));
//            m_cacheMC.insert(sim1[i].address(), {A.getValNoCache(sim1[i]), B.getValNoCache(sim2[i]), C.getValNoCache(sim1[i]), D.getValNoCache(sim2[i])});
        }
        m_pc1.setEvents(list1);
        m_pc1.prepareCache();
        m_pcMC1.setEvents(sim1);
        m_pcMC1.prepareCache();
        if (m_sameTag){
            m_pc2.setEvents(list2);
            m_pc2.prepareCache();
            m_pcMC2.setEvents(sim2);
            m_pcMC2.prepareCache();
        }
    }

    const real_t normFromCache() const {
        real_t n0=m_Anorm * m_Bnorm + m_Cnorm * m_Dnorm;
        real_t n = 0;
        real_t N=0;
        int id=0;
        for (auto& [addresses, values] : m_cacheMC){
            N++;
            auto a = values[0];
            auto b = values[1];
            auto c = values[2];
            auto d = values[3];
            auto f = m_pcMC1.getValCache(id);
            if (m_sameTag) f -= m_pcMC2.getValCache(id);
            id++;
            n += -2 * std::real(a * b * std::conj(c * d) * exp(complex_t(0,f)));
        }
        
        return n0 + n/N;
    } 
    const real_t LLFromCache() const {
        real_t n = normFromCache();
        real_t ll =0 ;
        int id=0;
        for (auto& [addresses, values] : m_cache){
            auto f = m_pc1.getValCache(id);
            if (m_sameTag) f -= m_pc2.getValCache(id);
            id++;      
            auto pdf = std::norm(values[0] * values[1] * exp(complex_t(0, f/2)) - values[2] * values[3] * exp(complex_t(0, -f/2)))/n;
            ll += log(pdf);
        }
        return -2*ll;
    }

    TNtuple * dumpVals(std::string tagName){
        
        TNtuple * tup = new TNtuple( (tagName + "_vals").c_str(), (tagName + "_vals").c_str(), "aR:aI:bR:bI:cR:cI:dR:dI:dd");
        for (int i=0; i < (*m_events1).size(); i++){
            if (m_debug) std::cout<<"\rat "<<i;
            auto v = getVals((*m_events1)[i], (*m_events2)[i]);
            auto a = v[0];
            auto b = v[1];
            auto c = v[2];
            auto d = v[3];
            real_t dd = std::imag(log(a * std::conj(c)/ (std::abs(a * std::conj(c)))  ));
            tup->Fill(std::real(a), std::imag(a),std::real(b), std::imag(b),std::real(c), std::imag(c),std::real(d), std::imag(d), dd);
        }
        return tup;
    }
    TNtuple * dumpValsMC(std::string tagName){
        
        TNtuple * tup = new TNtuple( (tagName + "_vals").c_str(), (tagName + "_vals").c_str(), "aR:aI:bR:bI:cR:cI:dR:dI:dd");
        for (int i=0; i < (*m_sim1).size(); i++){
            auto v = getVals((*m_sim1)[i], (*m_sim2)[i]);
            auto a = v[0];
            auto b = v[1];
            auto c = v[2];
            auto d = v[3];
            real_t dd = std::imag(log(a * std::conj(c)/ (std::abs(a * std::conj(c)))  ));
            tup->Fill(std::real(a), std::imag(a),std::real(b), std::imag(b),std::real(c), std::imag(c),std::real(d), std::imag(d), dd);
        }
        return tup;
    }
    void updateNorm(){
        m_norm = norm();
    }

    void prepareACBD(){
        m_ACstMC = {};
        m_BDstMC = {};
        for (size_t i=0;i<(*m_sim1).size();i++){
            m_ACstMC.push_back(m_A.getValNoCache((*m_sim1)[i]) * std::conj(m_C.getValNoCache((*m_sim1)[i])));
            m_BDstMC.push_back(m_B.getValNoCache((*m_sim2)[i]) * std::conj(m_D.getValNoCache((*m_sim2)[i])));
            
            m_sig_s01MC.push_back((*m_sim1)[i].s(0,1));
            m_sig_s02MC.push_back((*m_sim1)[i].s(0,2));


            if (m_sameTag){
                m_tag_s01MC.push_back((*m_sim2)[i].s(0,1));
                m_tag_s02MC.push_back((*m_sim2)[i].s(0,2));
            }

        }
    }

    complex_t getValForNorm(const int i) const;

    protected:
        double  m_norm  =    {0};
        MinuitParameterSet m_mps;
        std::string m_SFType;
        
        
        
        CoherentSum m_A ;
        CoherentSum  m_B;
        CoherentSum  m_C;
        CoherentSum  m_D;


        EventList *  m_events1 = {nullptr};
        EventList * m_sim1  = {nullptr};
        EventList * m_events2 = {nullptr};
        EventList * m_sim2  = {nullptr};
        
        EventType m_type1;
        EventType m_type2;   

        real_t m_normA;
        real_t m_normB;
        real_t m_normC;
        real_t m_normD;
        complex_t m_zAC; 
        complex_t m_zBD; 
        std::vector<std::vector<double> > m_poly;
        std::vector<FitFraction> outputFractionsA;
        std::vector<FitFraction> outputFractionsB; 
        std::vector<FitFraction> outputFractionsC;
        std::vector<FitFraction> outputFractionsD;
     
    
        
    
        bool m_debug;
        int m_debugFreq;
        bool m_pdebug;
        bool m_pNorm;
        bool m_updateNorms;
        int m_order;
        int m_orderMag;
        std::string m_polyType;
        std::string m_polyTypeMag;
        size_t m_prepareCalls = 0;
        real_t m_Anorm = 0;
        real_t m_Bnorm = 0;
        real_t m_Cnorm = 0;
        real_t m_Dnorm = 0;
        bool m_analyticNorm=false;
        bool m_sameTag;
        bool m_slowNorm;
        bool m_flat;
        PhaseCorrection m_pc1;
        PhaseCorrection m_pcMC1;
        PhaseCorrection m_pc2;
        PhaseCorrection m_pcMC2;
        complex_t m_BDMCSum;
        std::vector<complex_t> m_ACstMC = {};
        std::vector<complex_t> m_BDstMC = {};

        std::map<std::vector<real_t *> , std::vector<complex_t> > m_cache = {};
        std::map<std::vector<real_t *>, std::vector<complex_t> > m_cacheMC = {};

        std::vector<double> m_sig_s01MC;
        std::vector<double> m_tag_s01MC;
        std::vector<double> m_sig_s02MC;
        std::vector<double> m_tag_s02MC;

        bool m_calcZAtStart;

  };
}
#endif
