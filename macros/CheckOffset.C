void CheckOffset(TH1D * h){
  double x = 0;
  double y = 0;
  double cnt = 0;
  for(int i = 1; i<=h->GetNbinsX(); i++) {
    double ang = h->GetXaxis()->GetBinLowEdge(i);
    x+=h->GetBinContent(i)*cos(2*ang);
    y+=h->GetBinContent(i)*sin(2*ang);
    cnt+=h->GetBinContent(i);
  }
  cout<<"x: "<<x/cnt<<endl;
  cout<<"y: "<<y/cnt<<endl;
}
