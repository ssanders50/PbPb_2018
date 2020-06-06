void resubmit(){
  FILE * fin = fopen("copy.log","r");
  char buf[200];
  bool foundFailed = false;
  bool found100 = false;
  while(fgets(buf,200,fin)!=NULL) {
    string schk = buf;
    if(schk.find("failed")!=std::string::npos) {
      if(schk.find("status")!=std::string::npos){
	//cout<<schk;
	foundFailed=true;
      }
    }
    if(schk.find("finished")!=std::string::npos) {
      if(schk.find("100.0")!=std::string::npos){
	//cout<<schk;
	found100 = true;
      }
    }
    if((foundFailed||found100) && schk.find("Log file")!=std::string::npos) {
      int iloc = schk.find("crab");
      int ilocEnd = schk.find("/crab.log");
      string substr = schk.substr(iloc,ilocEnd-iloc+1);
      string sub;
      if(foundFailed) {
	sub = "crab resubmit --dir="+substr;
      } else if (found100) {
	sub = "rm -rf "+substr;
      }
      cout<<" ======================================== "<<endl;
      cout<<sub<<endl;
      system(sub.data());
      foundFailed=false;
      found100= false;
    }
  }
}
