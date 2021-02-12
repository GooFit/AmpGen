#include <TH1D.h>

void test_root()
{
  TH1D* h = new TH1D("h","", 100,-5,5);
  for( int i = 0 ; i != 1000; ++i ) h->Fill( gRandom->Gaus(0,1) );
  h->Draw();
  gPad->SaveAs("test_root.pdf");
}
