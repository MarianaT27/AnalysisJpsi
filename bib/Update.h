#ifndef Update
#define Update



void addLeptonIDscore(string name_file) { 

	//Electron variables
    Double_t electron_p, electron_theta, electron_phi;
    Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
    Double_t electron_sfpcal, electron_sfecin, electron_sfecout;

    //Positron variables
    Double_t positron_p, positron_theta, positron_phi;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;

    cout<<"Updating tree"<<endl;

	TString root_file = "/lustre19/expphy/volatile/clas12/mtenorio/Root/"+name_file+".root";
	TFile *f = new TFile(root_file,"update"); 
	TTree *tree = (TTree*)f->Get("analysis");
	TMVA::Reader *readerTMVA = new TMVA::Reader( "!Color:!Silent" );
	int model=9;
    // Create a set of variables and declare them to the reader
	Float_t P, Theta, Phi, PCAL,ECIN,ECOUT;

	Double_t score_pos_6,score_ele_6;
	Float_t m2PCAL=-1;
	Float_t m2ECIN=-1;
	Float_t m2ECOUT=-1;
	Float_t Nphe;


	//readerTMVA->AddVariable( "P",&P );
	//readerTMVA->AddVariable( "Theta",&Theta);
	//readerTMVA->AddVariable( "Phi",&Phi);
     //readerTMVA->AddVariable( "Nphe",&Nphe);
	readerTMVA->AddVariable( "SFPCAL",&PCAL);
	readerTMVA->AddVariable( "SFECIN",&ECIN);
	readerTMVA->AddVariable( "SFECOUT",&ECOUT );
	readerTMVA->AddVariable( "m2PCAL",&m2PCAL);
	readerTMVA->AddVariable( "m2ECIN",&m2ECIN);
	readerTMVA->AddVariable( "m2ECOUT",&m2ECOUT);

    //Book Methods
	TString weightfile_pos; 
	TString weightfile_ele;

	weightfile_ele= "/work/clas12/mtenorio/ML_weights_pass2/"+name_file+"neg/TMVAClassification_BDT_6.weights.xml";
	weightfile_pos= "/work/clas12/mtenorio/ML_weights_pass2/"+name_file+"pos/TMVAClassification_BDT_6.weights.xml";

	readerTMVA->BookMVA( "BDT pos method", weightfile_pos );
	readerTMVA->BookMVA( "BDT ele method", weightfile_ele );

	TBranch *score_p_6 = tree->Branch("score_pos_6",&score_pos_6,"score_pos_6/d");
	TBranch *score_e_6 = tree->Branch("score_ele_6",&score_ele_6,"score_ele_6/d");

	cout<<"Adding new branches"<<endl;

	tree->SetBranchAddress("electron_p",&electron_p);
	tree->SetBranchAddress("electron_theta",&electron_theta);
	tree->SetBranchAddress("electron_phi",&electron_phi);
	tree->SetBranchAddress("electron_m2ecin",&electron_m2ecin);
	tree->SetBranchAddress("electron_m2ecout",&electron_m2ecout);
	tree->SetBranchAddress("electron_m2pcal",&electron_m2pcal);
	tree->SetBranchAddress("electron_sfecin",&electron_sfecin);
	tree->SetBranchAddress("electron_sfpcal",&electron_sfpcal);
	tree->SetBranchAddress("electron_sfecout",&electron_sfecout);

	tree->SetBranchAddress("positron_p",&positron_p);
	tree->SetBranchAddress("positron_theta",&positron_theta);
	tree->SetBranchAddress("positron_phi",&positron_phi);
	tree->SetBranchAddress("positron_m2ecin",&positron_m2ecin);
	tree->SetBranchAddress("positron_m2ecout",&positron_m2ecout);
	tree->SetBranchAddress("positron_m2pcal",&positron_m2pcal);
	tree->SetBranchAddress("positron_sfecin",&positron_sfecin);
	tree->SetBranchAddress("positron_sfpcal",&positron_sfpcal);
	tree->SetBranchAddress("positron_sfecout",&positron_sfecout);


	Long64_t nentries = tree->GetEntries(); 
	cout<<"Starting loop"<<endl;
	for (Long64_t i=0;i<nentries;i++) {
		tree->GetEntry(i); 
		P=positron_p;
		Theta=positron_theta/57.2958;
		Phi=positron_phi/57.2958;
		PCAL=positron_sfpcal;
		ECIN=positron_sfecin;
		ECOUT=positron_sfecout;
		m2PCAL=positron_m2pcal;
		m2ECIN=positron_m2ecin;
		m2ECOUT=positron_m2ecout;
        score_pos_6=readerTMVA->EvaluateMVA("BDT pos method");

		P=electron_p;
		Theta=electron_theta/57.2958;
		Phi=electron_phi/57.2958;
		PCAL=electron_sfpcal;
		ECIN=electron_sfecin;
		ECOUT=electron_sfecout;
		m2PCAL=electron_m2pcal;
		m2ECIN=electron_m2ecin;
		m2ECOUT=electron_m2ecout;
        score_ele_6=readerTMVA->EvaluateMVA("BDT ele method");

		score_p_6->Fill(); 
		score_e_6->Fill();
	} 
	tree->Print(); 
	tree->Write(); 
	delete f; 
}

void addHelicity(string name_file, int version){
	 //Variables
    int run_hel,event_hel,hel_hel;

	cout<<"Creating tree"<<endl;

	TString root_file = "/lustre19/expphy/volatile/clas12/mtenorio/Root/"+name_file+"_hel.root";
    TFile *file = new TFile(root_file,"RECREATE");

	//TTree *tree = (TTree*)file->Get("analysis");
	TTree *helicity_tree = new TTree("analysis",root_file);

	helicity_tree->Branch("run_hel",&run_hel,"run_hel/i");
    helicity_tree->Branch("event_hel",&event_hel,"event_hel/i");
    helicity_tree->Branch("hel_hel",&hel_hel,"hel_hel/I");

	int max;
    if(version==-18)
        max=170;
     if(version==+18)
        max=180;
    if(version==-19)
        max=120;

	for(int fc =0; fc<max; fc++) {//Run 5032 to 5419 // 6616 6783

			char filename1[500];

			if(version==-18)
				sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/jpsitcs/jpsitcs_00%d.hipo",F18in_P2[fc]);
			else if(version==+18)
				sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/jpsitcs/jpsitcs_00%d.hipo",F18out_P2[fc]); 
			else if(version==-19)
				sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/jpsitcs/jpsitcs_00%d.hipo",S19_P2[fc]);


			//cout<<"---------ADD Hipo reader and file -------"<<endl;
			hipo::reader  reader;
			reader.open(filename1);

			//cout<<"---------ADD and read dictonary -------"<<endl;
			hipo::dictionary  factory;
			reader.readDictionary(factory);

			//cout<<"---------Show factory -------"<<endl;
			//factory.show();

			//cout<<"---------Call event -------"<<endl;
			hipo::event      event;

			hipo::bank EVENT(factory.getSchema("REC::Event"));
			hipo::bank HEADER(factory.getSchema("RUN::config"));
			hipo::bank PART(factory.getSchema("REC::Particle"));

			int counter = 0;
			while(reader.next()==true ){//Loops all events


				reader.read(event);

				event.getStructure(EVENT);
				event.getStructure(HEADER);
				event.getStructure(PART);


				int rn = 0;
				int en = 0;

				if(PART.getSize()<1) {
				continue;
				}

				if(HEADER.getRows()==1) {
					for(int i = 0; i < HEADER.getRows(); i++) {
						rn = HEADER.getInt("run",i);
						en = HEADER.getInt("event",i);
					}
				}
				
				auto hel = 0 ;
				for(int i = 0; i < EVENT.getRows(); i++) {
					hel= EVENT.getByte("helicity",i);
				}
			

				int number_of_electrons = 0;
				int number_of_positrons = 0;

				for(int i = 0; i < PART.getRows(); i++){
					int   pid = PART.getInt("pid",i);
					int    status = PART.getInt("status",i);
					//Status Between 2000 and 4000 Means Forward Detector
					if(abs(status)>=2000 && abs(status)<4000) {
						//PID 11 Means Electron
						if(pid==11) {
							number_of_electrons = number_of_electrons + 1;
						}
						//PID -11 Means Positron
						else if(pid==-11) {
							number_of_positrons = number_of_positrons + 1;
						}
					}
				}
		
				if(number_of_electrons==1 && number_of_positrons==1){			
					run_hel =rn;
					event_hel=en;
					hel_hel=hel;
					helicity_tree->Fill();
				}
			
			}//end while
			printf("Version= %d  run = %d\n",version, fc);
		}//End for "Runs"

	//tree->Print(); 
	file->Write(); 
}


#endif
