
#ifndef PCOHERENTSUM
#define PCOHERENTSUM


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
   pCoherentSum is the Quantum Correlated sum, designed as an input for the QcFitter program. 
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


class pCoherentSum {
    public:
      //Takes amplitudes - 
        pCoherentSum();
        pCoherentSum(const EventType& type1, const MinuitParameterSet& mps, 
        int gammaSign = 1, bool useXY = true, bool BConj = false );

        virtual ~pCoherentSum()=default;
        real_t operator()( const Event& event1) const { return prob(event1); }
        real_t prob(const Event& event1 ) const {
        double P = std::norm(getVal(event1))/m_norm;  
    
    if (m_debug){
        double A2 = std::norm(m_A.getVal(event1))/m_norm;

        double C2 = std::norm(m_C.getVal(event1))/m_norm;

        complex_t AC = m_A.getVal(event1) * std::conj(m_C.getVal(event1))/m_norm;

        auto i = Constant(0,1);
        complex_t eif = exp(i() * correction(event1));
        double inter = -2*std::real(AC);
        double corrected_inter = -2*std::real(AC * eif);       

    }
    return P;
        }

    double getC(int i, int j, std::string pref="")const {
        std::string key = "pCoherentSum::C"+pref+std::to_string(i)+std::to_string(j);
        double val=0;
        //if (m_debug) INFO(m_mps[key]->name()<<" = "<<m_mps[key]->mean());
        val = m_mps[key]->mean(); 
        return val;
    }
    double getdC2(int i, int j, std::string pref="")const {
        std::string key = "pCoherentSum::C"+pref+std::to_string(i)+std::to_string(j);
        double val=0;
        //if (m_debug) INFO(m_mps[key]->name()<<" = "<<m_mps[key]->mean());
        val = pow(m_mps[key]->err(), 2); 
        return val;
    }



    complex_t getSumFactor()const {

    if (m_useXY){
        if (m_gammaSign==-1){
            return complex_t(m_mps["pCoherentSum::x-"]->mean(), m_mps["pCoherentSum::y-"]->mean());
        }
        if (m_gammaSign==1){
            return complex_t(m_mps["pCoherentSum::x+"]->mean(), m_mps["pCoherentSum::y+"]->mean());
        }
    }
    else{
        return m_mps["pCoherentSum::rB"]->mean() * exp(complex_t(0, m_mps["pCoherentSum::deltaB"]->mean() + m_gammaSign * m_mps["pCoherentSum::gamma"]->mean()));
        

    }

        
    }
    
    CoherentSum getA(){
        return m_A;
    } 

    CoherentSum getC(){
        return m_C;
    } 

    const complex_t A(const Event& evt) const{
        return m_A.getVal(evt);
    }

    const complex_t C(const Event& evt) const{
        return m_C.getVal(evt);
    }
    real_t testnorm(){
        INFO("norm = "<<norm());
        return norm();
    }



    real_t prob_unnormalised(const Event& event1 ) const {return std::norm(getVal(event1));}
    void prepare();
    void reset(bool resetEvents);
    double slowNorm();
    void setEvents(EventList& list1);
    void setMC(EventList& list1);
    void setMC1(EventList& list1);


    void updateNorms(const std::vector<size_t>& iA, const std::vector<size_t>& iB,
        const std::vector<size_t>& iC, const std::vector<size_t>& iD);

    void debugNorm();
    void debug(const Event& event1) const; 
    complex_t getVal(const Event& event1) const;
    complex_t getValNoCache(const Event& event1) const;
    complex_t getValNoCache(const Event& event1, const size_t& offset) const;
    real_t size()const { return m_A.size() + m_C.size() ;}
    std::vector<std::vector<FitFraction> > fitFractions(const LinearErrorPropagator& linProp);
    real_t norm()  const;
    double probA(const Event& event){
        double prob = m_A.prob(event);
        return prob;
    }
    std::vector<complex_t> getVals(const Event& event1) const {
        complex_t A = m_A.getVal(event1);

        complex_t C = m_C.getVal(event1);

        complex_t ABCD = getVal(event1);///std::sqrt(m_norm);
        complex_t corr = correction(event1);
        std::vector<complex_t> vals = {A,C,ABCD, corr};
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
                real_t Cij = m_mps["pCoherentSum::C"+std::to_string(i) + std::to_string(j)]->mean();
                sum_i += Cij*fcn::legendre(zp,i) * fcn::legendre(zm,2*j+1);
               
                }
                corr += sum_i;
    }
    //corr = Constant(0,1) * corr;
    //complex_t val = exp(corr());
    if (m_pdebug) INFO("correction = "<<corr());
    return corr();

 
    }

    void updateNorm(){
        m_norm = norm();
    }

real_t norm_manual() const{
 if (m_debug) INFO("Getting the value for the normalised pdf");
  //complex_t sumFactor = getSumFactor(); 
  complex_t z(0,0);
  complex_t sumFactor = getSumFactor(); 
 for (size_t i=0; i < m_sim1.size(); i++){
    real_t f = m_pcMC1.calcCorrL((m_sim1)[i]);

    z += (m_A.getVal((m_sim1)[i]) * std::conj(m_C.getVal((m_sim1)[i]) * sumFactor ) * exp(complex_t(0,f)));
    z = z/(real_t)m_sim1.size(); 
 }

 return m_A.norm() + m_C.norm() + 2 * std::real(z);

}

real_t LL(){
    real_t _LL = 0;
    real_t n = norm();

    #pragma omp parallel for reduction( +: _LL )
    for (size_t i=0; i < m_events1.size(); i++){
        _LL += log(std::norm(getVal((m_events1)[i]))/n);
    }
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
                real_t Cij = m_mps["pCoherentSum::C"+std::to_string(i) + std::to_string(j)]->mean();
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

    void prepareCache( EventList& list1,   EventList& sim1)  {
        CoherentSum A(m_type1, m_mps);       
        
        A.transferParameters();
        A.setEvents(list1);
        A.setMC(sim1);
        A.prepare();

        auto Anorm = m_A.norm();
        INFO("m_ANorm = "<<A.norm());
        INFO("m_ANorm = "<<Anorm);
        CoherentSum C(m_type1.conj(true), m_mps);  C.transferParameters(); C.setEvents(list1); C.setMC(sim1); C.prepare(); m_Cnorm = m_C.norm();



        INFO("m_CNorm = "<<C.norm());


        for (size_t i=0; i<list1.size(); i++){
            complex_t a(A.getValNoCache(list1[i]));

            complex_t c(C.getValNoCache(list1[i]));

            m_cache.insert(std::pair<real_t* , std::vector<complex_t> >(list1[i].address(),  std::vector<complex_t> ({a, c})));
//            m_cache.insert(list1[i].address(), std::vector<complex_t>({A.getValNoCache(list1[i]), B.getValNoCache(list2[i]), C.getValNoCache(list1[i]), D.getValNoCache(list2[i])}));
//            m_cacheMC.insert(sim1[i].address(), {A.getValNoCache(sim1[i]), B.getValNoCache(sim2[i]), C.getValNoCache(sim1[i]), D.getValNoCache(sim2[i])});
        }
        for (size_t i=0; i<sim1.size();i++){
            m_Anorm += std::norm(A.getValNoCache(sim1[i]))/sim1.size();

            m_Cnorm += std::norm(C.getValNoCache(sim1[i]))/sim1.size();

        }
        INFO("a = "<<A.prob(list1[0]));
        INFO("m_ANorm = "<<m_Anorm);

        INFO("m_CNorm = "<<m_Cnorm);


        for (size_t i=0; i<sim1.size(); i++){
            complex_t a(A.getValNoCache(sim1[i]));

            complex_t c(C.getValNoCache(sim1[i]));


//            m_cacheMC.insert(std::pair<real_t*, std::vector<complex_t> >(sim1[i].address(),  std::vector<complex_t> ({a, b, c, d})));
            m_cacheMC.insert(std::pair<real_t* , std::vector<complex_t> >(sim1[i].address(),  std::vector<complex_t> ({a,  c })));
//            m_cache.insert(sim1[i].address(), std::vector<complex_t>({A.getValNoCache(sim1[i]), B.getValNoCache(sim2[i]), C.getValNoCache(sim1[i]), D.getValNoCache(sim2[i])}));
//            m_cacheMC.insert(sim1[i].address(), {A.getValNoCache(sim1[i]), B.getValNoCache(sim2[i]), C.getValNoCache(sim1[i]), D.getValNoCache(sim2[i])});
        }
        m_pc1.setEvents(list1);
        m_pc1.prepareCache();
        m_pcMC1.setEvents(sim1);
        m_pcMC1.prepareCache();
    }

    const real_t normFromCache() const {
        real_t n0=m_Anorm  + m_Cnorm ;
        real_t n = 0;
        real_t N=0;
        complex_t sf = getSumFactor();
        int id=0;
        for (auto& [addresses, values] : m_cacheMC){
            N++;
            auto a = values[0];
//            INFO("a = "<<a);

            auto c = values[1];

//            INFO("c = "<<c);

            auto f = m_pcMC1.getValCache(id);
            id++;
//            INFO("f = "<<f);

            n += 2  * std::real(a  * std::conj(c*sf) * exp(complex_t(0,f)));
        }
        
        return n0 + n/N;
    } 
    const real_t LLFromCache() const {
     //   INFO("Calculating norm from cache");
        real_t n = normFromCache();
        INFO("n = "<<n);
        real_t ll =0 ;
        complex_t sf = getSumFactor();
        int id =0;
        for (auto& [address, values] : m_cache){
            auto f = m_pc1.getValCache(id);
            id++;

            auto pdf = std::norm(values[0] * exp(complex_t(0, f/2)) + sf* values[1] * exp(complex_t(0, -f/2)))/n;
            ll += log(pdf);
        }
        return -2*ll;
    }


    TNtuple * dumpVals(std::string tagName){
        
        TNtuple * tup = new TNtuple( (tagName + "_vals").c_str(), (tagName + "_vals").c_str(), "aR:aI:cR:cI:dd");
        for (int i=0; i < m_events1.size(); i++){
            auto v = getVals(m_events1[i]);
            auto a = v[0];
            auto c = v[1];


            real_t dd = std::imag(log(a * std::conj(c)/ (std::abs(a * std::conj(c)))  ));
            tup->Fill(std::real(a), std::imag(a), std::real(c), std::imag(c), dd);
        }
        return tup;
    }
    TNtuple * dumpValsMC(std::string tagName){
        
        TNtuple * tup = new TNtuple( (tagName + "_vals").c_str(), (tagName + "_vals").c_str(), "aR:aI:cR:cI:dd");
        for (int i=0; i < m_sim1.size(); i++){
            auto v = getVals(m_sim1[i]);
            auto a = v[0];
            auto c = v[1];

            real_t dd = std::imag(log(a * std::conj(c)/ (std::abs(a * std::conj(c)))  ));
            tup->Fill(std::real(a), std::imag(a), std::real(c), std::imag(c),dd);
        }
        return tup;
    }
   
    protected:
        double  m_norm  =    {0};
        MinuitParameterSet m_mps;

        CoherentSum  m_A;

        CoherentSum  m_C;



        EventList m_events1;
        EventList m_sim1 ;


        
        EventType m_type1;


        real_t m_normA;

        real_t m_normC;

      
        std::vector<std::vector<double> > m_poly;
        std::vector<FitFraction> outputFractionsA;

        std::vector<FitFraction> outputFractionsC;

     
    
        
    
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

        real_t m_Cnorm = 0;

        bool m_sameTag;
        bool m_slowNorm;
        bool m_flat;
        PhaseCorrection m_pc1;
        PhaseCorrection m_pcMC1;

        int m_gammaSign;
        bool m_useXY;
        bool m_BConj;


        std::map<real_t * , std::vector<complex_t> > m_cache = {};
        std::map<real_t *, std::vector<complex_t> > m_cacheMC = {};



  };
}
#endif
