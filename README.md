# PbPb_2018
Event plane analysis code for CMS 2018 PbPb run

cmsrel CMSSW_10_3_1_patch3   
cd CMSSW_10_3_1_patch3/src  
cmsenv  
git cms-init  
git remote add CmsHI git@github.com:CmsHI/cmssw.git  
git fetch CmsHI  
git checkout CmsHI/forest_CMSSW_10_3_1 -b forest_CMSSW_10_3_1  
git cms-addpkg RecoVertex/PrimaryVertexProducer  
git cms-addpkg HeavyIonsAnalysis  

