// grp2.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#define ld long double

/* Globals to be used multiple times in the code */
// total recorded files available : 30
int numTrainfiles = 20, numTestfiles = 0;
string ROLL = "160101050";
string ORIGINAL_ROLL = "160101050";
string ROLL1 = "160101050";
string ROLL2 = "160101001";
string DEFAULT_ROLL = "NEW_USER";
string Recordings_Dir = ROLL+"_Recordings/"; // last string character must be '/'
string genUniversefile = "genUniverse.txt";
string genCodeBookfile = "genCodeBook.txt";
string WORDS[] = {"0","1","2","3","4","5","6","7","8","9","BEGIN","NO","PAUSE","QUIT","RESTART","RESUME","YES"};
string CUR_WORD;
int TRAINABLE_MODULE_REC = 10;
string rollPlayer1=ROLL, rollPlayer2=ROLL;
struct stat info;
//////////////////////////////////////////////////////////// Compute_Ci //////////////////////////////////////////////////////////
const ld PI = 4*atan(1.0);
const int p = 12;   // "p+1" is the size of Ri, Ai, & Ci vectors
const int m = 320;  // batch frame size
const ld tokuraDistWeight[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
const int norm_max = 5000;   // normalising factor
const ld factor = 0.001;   // fractional factor value of the max ste value for selecting the required region frames and excluding rest only-noise-containing frames
const int SLIDE = 80;
/*
modified data : sample array of amplitude values which is modified later...
original data : original sample array of amplitude values
ste : vector to store the short term energy values per frame
max_ste : vector to store max short term energy value (max_ste[0]) and corresponding index (max_ste[1]) over all the frames
R_file, A_file, C_file : stores the Ri, Ai & Ci values resp. for the current file
*/

vector<ld> modified_data, original_data, ste, max_ste(2), max_min_amp, R_file, A_file, C_file;
vector<vector<ld>> Universe_Ci; // stores the Ci vectors for ALL the frames (inside marker regions) of all(30) recorded files 
vector<int> markers; // specifies the valid region of samples in a file without noise

// to find max and min amplitude value among all samples of a file
// max_amp --> indx(0); max_amp_indx --> indx(1); min_amp --> indx(2); min_amp_indx --> indx(3); 
void find_max_min(vector<ld> arr, int n){
	max_min_amp[0] = (ld)(-1.0);
	max_min_amp[2] = (ld)FLT_MAX;

	for (int i = 0; i < n; ++i)
	{
		if(max_min_amp[0] < arr[i]){
			max_min_amp[0] = arr[i];
			max_min_amp[1] = i;
		}
		if(max_min_amp[2] > arr[i]){
			max_min_amp[2] = arr[i];
			max_min_amp[3] = i;
		}
	} 

	//cout<<"Max_Amp : "<<max_min_amp[0]<<endl;
	//cout<<"Min_Amp : "<<max_min_amp[2]<<endl;
	return ;
}

// to normalise all amplitude values so that they are within [-norm_max, norm_max]
void normalise(vector<ld> &arr, int n, ld max_amp, ld min_amp){
	/*
	ld a = (ld)(-1.0)*(ld)(norm_max);
	ld b = (ld)(norm_max);
	cout<<a<<endl;
	cout<<max_amp<<endl;
	cout<<min_amp<<endl;
	for (int i = 0; i < n; ++i)
	{
		arr[i] = ((ld)((arr[i]-min_amp) * (b-a)) / (ld)(max_amp-min_amp)) + a;
	}  
	*/
	ld maxValue = max(abs(max_amp), abs(min_amp));
	for (int i = 0; i < n; ++i)
	{
		arr[i] = (ld)(arr[i] * (ld)(norm_max)) / maxValue;
	}

	return;
}

// to compute the ste values per frame alongwith max ste value and store in corresponding vector 
void compute_ste(vector<ld> arr, int n){
	max_ste[0] = -1;
	max_ste[1] = 0;
	ste.clear();

	for(int i=0; i<=n-m; i+=m){
		ld temp_ste = 0;

		// finding the ste and zcr values for the current frame
		for(int j=i; j<i+m; j++){
			temp_ste += (ld)((arr[j] * arr[j]) / (ld)(m) ) ;
		}

		// finding the max ste value
		if(max_ste[0] < temp_ste){
			max_ste[0] = temp_ste;
			max_ste[1]= i/m;
		}

		// storing values (per frame) in corresponding vector
		ste.push_back(temp_ste);
	}

	if(max_ste[0] == -1) cout<<endl<<endl<<"ERROR : Max_STE cannot be -1 !!"<<endl<<endl;

	return;
}

// finding markers by traversing either way from the frame with max ste value till the current frame ste value drops to the "factor" amount of the max_ste value 
void find_markers(vector<ld> arr,int n){
	ld comp_ste = factor * max_ste[0] ;

	// left_marker (starting from leftmost index of sample inside the concerned region)
	int lm = (int)max_ste[1] - 1;
	while(lm >= 0  &&  ste[lm] >= comp_ste) lm--;
	if(lm<0) lm=0;

	// right_marker (ending at right index of sample just outside the concerned region)
	int rm = (int)max_ste[1] + 1;
	while((rm < ste.size()) && ste[rm] >= comp_ste) rm++;
	if(rm == ste.size()) rm--;

	markers[0] = lm*m ;
	markers[1] = rm*m + m ;
	
	//cout<<"Markers : "<<markers[0]<<" "<<markers[1]<<"  "<<n<<endl;

	return;
}

// find avg noise in left_over portion outside marker regions and subtract that from every sample for computing dc shift 
void dc_shift(vector<ld> &arr, int n){
	// calculate average amplitude of noise
	ld avg_amp = 0;
	ld total_noise_data = n - (markers[1] - markers[0]);

	// iterating over left silent region
	for(int i=0; i< markers[0]; i++){
		avg_amp +=  (ld)(arr[i] / (total_noise_data)) ;
	}

	// iterating over right silent region
	for(int i=markers[1]; i<n; i++){
		avg_amp +=  (ld)(arr[i] / (total_noise_data)) ;
	}

	//cout<<"avg_noise_amp : "<<avg_amp<<endl;

	// subtract average amplitude of noise from the whole data
	// dc shift
	for(int i=0; i<n; i++){
		arr[i] -= avg_amp;
	}

	return;
}

// prints the given vector
void printarr(vector<ld> arr, int n){
	for(int i=0; i<n; i++){
		cout<<arr[i]<<" ";
	}
	cout<<endl<<endl<<endl;
}

// Finds the Ri vector for the given file
void find_ri(vector<ld> &arr, int n, vector<ld> &r){
	ld sum;
	for(int k=0; k<=p; k++){
		sum=0;
		for(int i=0; i<=n-1-k; i++){
			sum +=  (ld)(arr[i] * arr[i+k]);
		}
		r[k] = sum;
	}
}

// Finds the Ai vector for the given file
void find_ai(vector<ld> &a, vector<ld> r){
	vector<ld> e(p+1), k(p+1);
	vector<vector<ld>> v(p+1,vector<ld>(p+1));
	e[0] = r[0];
	k[1] = r[1]/r[0] ;

	for(int i=1; i<=p; i++){

		if(i!=1){
			k[i] = r[i]/e[i-1];
			for(int j=1; j<=i-1; j++){
				k[i] -= (r[i-j] * v[j][i-1]) / e[i-1];
			}

			for(int j=1; j<=i-1; j++){
				v[j][i] = v[j][i-1] - (k[i] * v[i-j][i-1]);
			}
		}
		v[i][i] = k[i];
		e[i] = e[i-1] * (1 - k[i]*k[i]) ;
	}

	for(int i=1; i<=p; i++){
		a[i] = v[i][p];
	}

}

// Finds the Ci vector for the given file
void find_ci(vector<ld> &c, vector<ld> a){
	c[1] = a[1];
	for(int i=2; i<=p; i++){
		c[i] = a[i];
		for(int j=1; j<=i-1; j++){
			c[i] += (ld)( ((ld)j * c[j] * a[i-j]) / (ld)(i));
		}
	}
}

// applies the hamming function for the given sample
void find_h(vector<ld> &arr, vector<ld> &w, int n){
	for(int i=0; i<n; i++){
		w[i] = 0.54 - (0.46 * cos( (2*PI*i) / (ld)(n-1)));
		arr[i] = arr[i]*w[i];
	}
}

// applies the Liftering Sine function for the given sample on the Ci values
void LiftererSineWindow(vector<ld> &C_arr){
	ld Wt;
	for(int j=0; j<(int)C_arr.size(); j++){
		Wt = 6 * sin(j * (PI / p))  +  1;
		C_arr[j] = C_arr[j] * Wt;
	}
}

// computes the Ri, Ai & Ci vectors for the file by taking the frames inside the marker region with a sliding window of length 80...
void compute_Ci_file(int n){
	    // for (int i=0; i<=n-m; i += SLIDE){
		for (int i=markers[0]; i<=markers[1]-m; i += SLIDE){
			int frame_id = i/SLIDE ;
			vector<ld> frame_data;

			for(int j = i; j<i+m; j++){
				frame_data.push_back(modified_data[j]);
			}

			vector<ld> R_frame(p+1), A_frame(p+1), C_frame(p+1);

			find_ri(frame_data,frame_data.size(),R_frame);
			R_file.insert(R_file.end(),R_frame.begin(),R_frame.end());

			find_ai(A_frame,R_frame);
			A_file.insert(A_file.end(),A_frame.begin(),A_frame.end());

			find_ci(C_frame,A_frame);
			LiftererSineWindow(C_frame);
			C_file.insert(C_file.end(),C_frame.begin(),C_frame.end());
		}
		return;
}

// dumps the Ci file in the appropriate folder
void dump_Ci_file(string fileName){
	string outputfileName;
	ofstream fout;
	outputfileName = "CiValues/" + fileName ;
	//cout<<outputfileName<<endl;
	fout.open(outputfileName);
	vector<ld> tmp;
	for(int id=0; id<(int)C_file.size(); id++){
		if(id%13){
			tmp.push_back(C_file[id]);
			fout<<C_file[id]<<" ";
		}
		else if(id!=0){
			Universe_Ci.push_back(tmp);
			tmp.clear();
			fout<<endl;
		}
	}
	// appending the current Ci vector to Universe
	Universe_Ci.push_back(tmp);
	fout.close();
}

// dumps the Trimmed Recording file of samples in the appropriate folder
void dumpTrimmedRec(string fileName){
	string outputfileName;
	ofstream fout;
	outputfileName = "TrimmedRecordings/" + fileName ;
	// cout<<outputfileName<<endl;
	fout.open(outputfileName);
	for(int id=markers[0]; id<markers[1]; id++){
		fout<<modified_data[id]<<endl;
	}
	fout.close();
}

// processes the given sample file, by applying dc shift, normalisation and finally computing the Ci values for all the frames inside the marker region
void processData(string fileNameFull, string fileName){
	ifstream file(fileNameFull);
	cout<<"fileNameFull : "<<fileNameFull<<endl;
	string str;
	int n=0;

	modified_data.clear();
	original_data.clear();

	//Read the Raw Data
	while(getline(file,str)){
		//cout<<"str : "<<str<<endl;
		if(str.length()==0) continue;

		if((str[0] >='0' && str[0] <= '9') || str[0]=='-'){
			n++;
			modified_data.push_back(stold(str));
			// cout<<"data : "<<modified_data.back()<<endl;
			original_data.push_back(stold(str));
		}
	}
	file.close();
	//cout<<"Data Size : "<<n<<endl;

	max_min_amp = vector<ld> (4);
	markers = vector<int> (2);

	find_max_min(modified_data,modified_data.size());
	normalise(modified_data,modified_data.size(),max_min_amp[0],max_min_amp[2]);
	ste.clear();
	compute_ste(modified_data,modified_data.size());
	find_markers(modified_data,modified_data.size());

	modified_data.clear();
	modified_data = original_data;

	dc_shift(modified_data, modified_data.size());
	find_max_min(modified_data,modified_data.size());
	normalise(modified_data,modified_data.size(),max_min_amp[0],max_min_amp[2]);
	dumpTrimmedRec(fileName);
	
	ste.clear();
	compute_ste(modified_data,modified_data.size());

	R_file.clear();
	A_file.clear();
	C_file.clear();
	compute_Ci_file(modified_data.size());
	dump_Ci_file(fileName);

	return;
}

// dumps the Universe Ci Vector to the file : genUniversefile
void dumpUniverseCi(){
	ofstream fout(genUniversefile);
	for(int i=0; i<Universe_Ci.size(); i++){
		for(int j=0; j<Universe_Ci[0].size(); j++){
			fout<<Universe_Ci[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
}

// processes the training & testing files, saves Ci vectors (& Trimmed Recordings) for each file in appropriate folder under "CiValues" folder (& "TrimmedRecordings" folder resp.) , and also dumps the universe Ci vector 
void Compute_Ci(int numfiles){

	string fileName, folderName, dir;

	dir = "CiValues";
	wstring wdir(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	dir = "TrimmedRecordings";
	wdir.clear();
	wdir.assign(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	string newname = CUR_WORD;
	vector<string> newfile;
	for(int i=1; i<=numfiles; i++){
		stringstream ss;
		ss << i;
		newfile.push_back(ss.str());
	}
	
	for(int k=0; k<1; k++){
		folderName = ROLL + "_" + newname ;
		dir = "CiValues/" + folderName;
		wstring wdir(dir.begin(), dir.end());
		if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
		else cout<<"False "<<endl;

		dir = "TrimmedRecordings/" + folderName;
		wdir.clear();
		wdir.assign(dir.begin(), dir.end());
		if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
		else cout<<"False "<<endl;

		for(int fid=0; fid<numfiles; fid++){
			fileName = folderName + "/" + ROLL + "_"+newname+"_"+newfile[fid]+".txt";

			R_file.clear();
			A_file.clear(); 
			C_file.clear();
			processData(Recordings_Dir + fileName, fileName);
		}
	}

	dumpUniverseCi();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////// genCodeBook /////////////////////////////////////////////////////////////
const int K = 32;  // final size of CodeBook to be generated
const ld deltaDiff = 0.001;  // value of delta difference b/w avg. distortion values for halting the K-Means algo.
ld newAvgDist=0, prevAvgDist = 1.0;  //  stores new(current) and prev avg. distortion relative to that iteration
ld eps = 0.05;  // default splitting parameter (currently not used)

// to store Universe vectors with their label of codebook vector they are associated to
typedef struct node{
	vector<ld> CiValues;
	int codeBookIndex;
} sn; 

vector<sn> universe_vectors;   // stores all the universe Ci vectors along with their labels
vector<vector<ld>> codebook_vectors;  // stores all the codebook vectors
sn temp_univ_vector;  // used for storing details of a single universe Ci vector
vector<ld> temp_CB_vector;  // used for storing details of a single codeBook vector


// assigns the label of codebook vector having minm Tokhura distance from the particular Universe vector of the given index (indx).
void calcMinDistance(int indx){
	ld min=FLT_MAX, cur=0.0, tmp;
	for(int i=0; i<codebook_vectors.size(); i++){
		cur=0.0;
		for(int j=0; j<codebook_vectors[0].size(); j++){
			tmp = (universe_vectors[indx].CiValues[j] - codebook_vectors[i][j]) ;
			cur +=  tokuraDistWeight[j]*tmp*tmp;
		}

		if(cur<min){
			min = cur;
			universe_vectors[indx].codeBookIndex = i;
		}
	}
	// cout<<min<<endl;
}

// assigns the codeBook label to each Universe vector
void assignCB(){
	for(int i=0; i<universe_vectors.size(); i++){
		calcMinDistance(i);
	}
}

// calculates the Avg. distortion over all the Universe vectors
void calcAvgDistortion(){
	ld tmp;
	prevAvgDist = newAvgDist;
	newAvgDist = 0.0;
	
	for(int i=0; i<universe_vectors.size(); i++){
		for(int j=0; j<universe_vectors[0].CiValues.size(); j++){
			tmp = (universe_vectors[i].CiValues[j] - codebook_vectors[universe_vectors[i].codeBookIndex][j]) ;
			newAvgDist +=  tokuraDistWeight[j]*tmp*tmp;
		}
	}

	newAvgDist = newAvgDist / (ld)(universe_vectors.size()) ;
}

// updates the CodeBook by assigning the new codebook vectors to be the mean of all the universe vectors assocaited to  that cluster
void updateCB(){
	// "counter" vector stores the count of universe vectors in each cluster
	vector<int> counter(codebook_vectors.size(), 0);
	
	// "sumXi" stores the sum of Ci values of universe vectors which are labeled to each codebook vector, for each vector in codeBook
	vector<vector<ld>> sumXi(codebook_vectors.size(), vector<ld>(codebook_vectors[0].size(), 0.0));
	
	for(int i=0; i<universe_vectors.size(); i++){
		for(int j=0; j<codebook_vectors[0].size(); j++){
			sumXi[universe_vectors[i].codeBookIndex][j] += universe_vectors[i].CiValues[j] ; 
		}
		counter[universe_vectors[i].codeBookIndex]++;
	}

	for(int i=0; i<codebook_vectors.size(); i++){
		// escaping the Empty Cell condition
		if(counter[i] != 0){
			for(int j=0; j<codebook_vectors[0].size(); j++){
				codebook_vectors[i][j] = sumXi[i][j] / (ld)(counter[i])  ;
			}
		}
	}
}

// implementing K-Means Algorithm
void kMeans(){
	int cnt=0;
	while(abs(prevAvgDist - newAvgDist)  > deltaDiff ){
		cnt++;
		assignCB();
		calcAvgDistortion();
		updateCB();
		//cout<<"NEW Avg. Distort. : "<<newAvgDist<<endl;
		//cout<<"PREV Avg. Distort. : "<<prevAvgDist<<endl;
		//cout<<endl<<endl;
	}
	cout<<"Iterations Count : "<<cnt<<endl<<endl;
}

// Calculates the splitting epsilon value for every centroid(codebook vector) 
vector<ld> calcSplit(vector<ld> &mean, int label, vector<sn> &universe){
	int n = universe.size(), cols=universe[0].CiValues.size();
	ld sum_dev;
	vector<ld> split_epsilon(cols);
	for (int j = 0; j < cols; ++j){
		sum_dev=0;
		for (int i = 0; i < n; ++i){
			if(universe[i].codeBookIndex == label)
				sum_dev += pow(universe[i].CiValues[j] - mean[j],2);
		}		
		//split_epsilon[j] = (ld)(sqrt(sum_dev/n)) / (ld)(1000.0);
		split_epsilon[j] = (ld)(sum_dev/n) / (ld)(1000.0);
	}
	return split_epsilon;
}

// generates the CodeBook using the Universe Ci vectors in the given file with name as : "fileName" 
void genCodeBook(string fileName){
	string word;
    ifstream file(fileName);

	universe_vectors.clear();
	temp_univ_vector.CiValues.clear();
	int cnt=0;  // takes 12 Ci values from each line for a Universe vector, gets reset to 0 after reaching 12

	while(file >> word){ 
		temp_univ_vector.CiValues.push_back(stod(word));
		cnt++;
		if(cnt==12){
			temp_univ_vector.codeBookIndex = 0;
			universe_vectors.push_back(temp_univ_vector);
			temp_univ_vector.CiValues.clear();
			cnt=0;
		}
	}
	file.close(); 

	cout<<"No. of vectors in Universe : "<<universe_vectors.size()<<endl<<endl;

	// pushing a dummy value as the first codeBook vector
	codebook_vectors.push_back(universe_vectors[0].CiValues);
	
	// assigns the mean of all the universe vectors as the codeBook vector 
	updateCB();

	// calculates initial avg. distortion value
	calcAvgDistortion();

	int m=1;  //  size of CodeBook

	//cout<<"NEW Avg. Distort. : "<<newAvgDist<<endl;
	//cout<<"PREV Avg. Distort. : "<<prevAvgDist<<endl<<endl<<endl;

	vector<ld> split_epsilon;

	while(m<K){
		cout<<"Size of Updated CodeBook : "<<m*2<<"   Final Size : "<<K<<endl;
		
		// doubling the codeBook 
		for(int i=0; i<m; i++){
			split_epsilon = calcSplit(codebook_vectors[i],i,universe_vectors);
			temp_CB_vector.clear();
			for(int j=0; j<codebook_vectors[0].size(); j++){
				//temp_CB_vector.push_back(codebook_vectors[i][j] * (1 + eps))  ;
				//temp_CB_vector.push_back(codebook_vectors[i][j] * (1 + split_epsilon[j]))  ;
				temp_CB_vector.push_back(codebook_vectors[i][j] + split_epsilon[j])  ;
			}

			for(int j=0; j<codebook_vectors[0].size(); j++){
				//codebook_vectors[i][j] = codebook_vectors[i][j] * (1 - eps)  ;
				//codebook_vectors[i][j] = codebook_vectors[i][j] * (1 - split_epsilon[j])  ;
				codebook_vectors[i][j] = codebook_vectors[i][j] - split_epsilon[j]  ;
			}

			codebook_vectors.push_back(temp_CB_vector);
		}

		kMeans();

		// doubling size for use in next iteration 
		m *= 2;

		// initialising prev. avg. dist. for the next doubled codeBook
		prevAvgDist=FLT_MAX;
	}

	ofstream fout(genCodeBookfile);
	for(int i=0; i<codebook_vectors.size(); i++){
		for(int j=0; j<codebook_vectors[0].size(); j++){
			fout<<codebook_vectors[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
	
	cout<<endl<<"Final Average Distortion : "<<newAvgDist<<endl<<endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////// GenerateObsvSeq //////////////////////////////////////////////////////////
vector<int> Obs_file; // for storing the Observation seq. for the whole file
vector<vector<ld>> codebook; // for storing the (32-sized) CodeBook Ci vectors

// finds and returns the Obs. seq. no. for the current frame
int getObsNo(vector<ld> Ci_frame){
	ld minDist = FLT_MAX, temp = 0;
	int minDistIndx = -2;
	for(int i=0; i<(int)codebook.size(); i++){
		temp = 0;
		for(int j=0; j<(int)codebook[0].size(); j++){
			temp += (tokuraDistWeight[j] * pow((Ci_frame[j] - codebook[i][j]),2)) ;
		}
		if(temp < minDist) minDist=temp, minDistIndx=i;
	}

	return minDistIndx+1;
}

// Generates the Obs. seq for the current file and populates the "Obs_file" vector
void genObsSeq(string fileName){
	string word;
	ifstream fin(fileName);
	vector<ld> Ci_frame;
	int cnt=0;
	while(fin>>word){
		cnt++;
		Ci_frame.push_back(stold(word));
		cnt = cnt % 12;
		if(cnt==0){
			Obs_file.push_back(getObsNo(Ci_frame));
			Ci_frame.clear();
		}
	}
	fin.close();
}

// Generates the Obs. seq for all the pre-recorded files and dumps each one of them in the appropriate folder under "ObsvSeq" folder
void GenerateObsvSeq(int numfiles, string codeBookFileName, bool rec_flg, string recfileName){
	string fileName, folderName, word, dir;
	dir = "ObsvSeq";
	wstring wdir(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	ifstream fin(codeBookFileName);
	codebook.clear();
	vector<ld> tmp;
	int cnt=0;

	while(fin>>word){
		cnt++;
		tmp.push_back(stold(word));
		cnt = cnt % 12;
		if(cnt==0){
			codebook.push_back(tmp);
			tmp.clear();
		}
	}
	fin.close();

	if(rec_flg==0){
		string newname = CUR_WORD;
		vector<string> newfile;
		for(int i=1; i<=numfiles; i++){
			stringstream ss;
			ss << i;
			newfile.push_back(ss.str());
		}
	
		for(int k=0; k<1; k++){
			folderName = ROLL + "_" + newname ;
			dir = "ObsvSeq/" + folderName;
			wstring wdir(dir.begin(), dir.end());
			if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
			else cout<<"False "<<endl;

			for(int fid=0; fid<numfiles; fid++){
				fileName = folderName +"/" + ROLL + "_" +newname+"_"+newfile[fid]+".txt";
				Obs_file.clear();
				genObsSeq("CiValues/" + fileName);

				ofstream fout("ObsvSeq/" + fileName);
				cout<<"ObsvSeq/"<<fileName<<endl;
				for(int num=0; num<(int)Obs_file.size(); num++){
					fout<<Obs_file[num]<<endl;
				}
				fout.close();
			}
		}
	}
	else{
		Obs_file.clear();
		genObsSeq("CiValues/" + recfileName);
		ofstream fout("ObsvSeq/" + recfileName);
		cout<<"ObsvSeq/"<<recfileName<<endl;
		for(int num=0; num<(int)Obs_file.size(); num++){
			fout<<Obs_file[num]<<endl;
		}
		fout.close();
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////// GenerateModels ///////////////////////////////////////////////////////////
int N=5;  // No. of States
int M=32; // No. of Different Observations (size of  codebook)
int T; // Total Obs. Time (from t=1 to t=T)
int TOTAL_ITER = 1; // max total iterations after which we want the output as the Avg. Model 
ld ProbOgivenLambda; // P(O/lambda)
ld Pstar; // max defining Probability for a model (for a particular iteration)
ld prevPstar; // defining Probability for a model (for the previous iteration)
ld BestPstar; // best defining Probability of a model (for a particular digit)
vector<ld> Pi(N+1); 
vector<vector<ld>> A(N+1,vector<ld>(N+1));
vector<vector<ld>> B(N+1,vector<ld>(M+1));
vector<ld> newPi(N+1);
vector<vector<ld>> newA(N+1,vector<ld>(N+1));
vector<vector<ld>> newB(N+1,vector<ld>(M+1));
vector<ld> AvgPi; 
vector<vector<ld>> AvgA;
vector<vector<ld>> AvgB;
vector<ld> BestPi; 
vector<vector<ld>> BestA;
vector<vector<ld>> BestB;
vector<int> O;  // current (selected) Obs. sequence
vector<int> Q;  // current State Sequence for the given Obs seq.
vector<vector<ld>> Alpha;
vector<vector<ld>> Beta;
vector<vector<ld>> Gamma;
vector<vector<ld>> Delta;
vector<vector<int>> Psi;
vector<vector<vector<ld>>> Xi;

// initialises (1-indexed --> ) A, B & Pi matrices  
void initialiseModel(string iterNo, string k){
	string word, fileNameA, fileNameB, fileNamePi;
	if(iterNo=="0"){
		fileNameA = "matrixA.txt";
		fileNameB = "matrixB.txt";
		fileNamePi = "matrixPi.txt";
	}
	else{
		fileNameA = "Models/" + ORIGINAL_ROLL + "_" + k + "_matrixA.txt";
		fileNameB = "Models/" + ORIGINAL_ROLL + "_" + k + "_matrixB.txt";
		fileNamePi = "Models/" + ORIGINAL_ROLL + "_" + k + "_matrixPi.txt";
	}


	ifstream fin;
	fin.open(fileNameA);
	//cout<<"Original A : "<<endl;
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fin >> word;
			A[i][j] = stold(word);
			//cout<<A[i][j]<<" ";
		}
		//cout<<endl;
	}
	fin.close();

	fin.open(fileNameB);
	//cout<<endl<<"Original B : "<<endl;        
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fin >> word;
			B[i][j] = stold(word);
			//cout<<B[i][j]<<" ";
		}
		//cout<<endl;
	}
	fin.close();

	fin.open(fileNamePi);
	//cout<<endl<<"Original Pi : "<<endl;
	for(int i=1; i<=N; i++){
		fin >> word;
		Pi[i] = stold(word);
		//cout<<Pi[i]<<" ";
	}
	fin.close();

	//cout<<endl<<endl<<endl;
}


// Solution to Problem-1 : Alpha, Beta
void prob1(){

	/******** Alpha ********/
	/* initialization for t=1 */
	//cout<<endl<<"Alpha : "<<endl;
	for(int i=1; i<=N; i++){
		Alpha[1][i] = Pi[i] * B[i][O[1]] ;
		//cout<<Alpha[1][i]<<" ";
	}
	//cout<<endl;

	/* induction for t>1 */
	for(int t=1; t<=T-1; t++){
		for(int j=1; j<=N; j++){
			Alpha[t+1][j] = 0 ;
			for(int i=1; i<=N; i++) 
				Alpha[t+1][j] += (Alpha[t][i] * A[i][j]) ;
			Alpha[t+1][j] *= B[j][O[t+1]] ;
			//cout<<Alpha[t+1][j]<<" ";
		}
		//cout<<endl;
		//if(Alpha[t+1][1] == 0) cout<<t+1<<" ZIROOOO "<<endl;
	}

	/* termination */ 
	ProbOgivenLambda = 0;
	for(int i=1; i<=N; i++) 
		ProbOgivenLambda += Alpha[T][i];

	//cout<<endl<<"ProbOgivenLambda : "<<ProbOgivenLambda<<endl;


	/******** Beta ********/
	/* initialization for all t (& hence t=T) done already while initializing Beta to be 1 */

	/* induction for t=T-1 to t=1 */
	for(int t=T-1; t>=1; t--){
		for(int j=1; j<=N; j++){
			Beta[t][j] = 0 ;
			for(int i=1; i<=N; i++) 
				Beta[t][j] += (Beta[t+1][i] * A[j][i] * B[i][O[t+1]]) ;
		}
	}

	
	//cout<<endl<<"Beta : "<<endl;
	
	for(int t=1; t<=T; t++){
		for(int j=1; j<=N; j++){
			//cout<<Beta[t][j]<<" ";
		}
		//cout<<endl;
		//if(Beta[t][1] == 0) cout<<t<<" ZIROOOO "<<endl;
	}
	
}


// Solution to Problem-2 : Gamma, Delta, Psi, Q
void prob2(){
	prevPstar = Pstar;
	/******** Gamma ********/
	//cout<<endl<<endl<<"Gamma : "<<endl;
	for(int t=1; t<=T; t++){
		ld sum = 0;
		for(int i=1; i<=N; i++){
			Gamma[t][i] =  (Alpha[t][i] * Beta[t][i]) ;
			sum += Gamma[t][i];
		}

		for(int i=1; i<=N; i++){
			if(sum!=0) Gamma[t][i] /=  sum ;
			//cout<<Gamma[t][i]<<" ";
		}
		//cout<<endl;
	}



	/******** Delta, Psi ********/
	/* initialization for t=1 */
	for(int i=1; i<=N; i++) 
		Delta[1][i] = Pi[i] * B[i][O[1]] ,    Psi[1][i] = 0 ;

	ld maxval = INT_MIN;
	int maxarg = INT_MIN;

	/* induction for t>1 */
	for(int t=2; t<=T; t++){
		for(int j=1; j<=N; j++){
			maxval = INT_MIN;
			maxarg = INT_MIN;
			Delta[t][j] = B[j][O[t]] ;
			for(int i=1; i<=N; i++) 
				if( (A[i][j] * Delta[t-1][i]) > maxval ) 
					maxval = (A[i][j] * Delta[t-1][i])  ,   maxarg = i;
			Delta[t][j] *= maxval ;
			Psi[t][j] = maxarg;
		}
	}

	/*
	cout<<endl<<endl<<"Delta : "<<endl;
	for(int t=1; t<=T; t++){
		for(int j=1; j<=N; j++){
			cout<<Delta[t][j]<<" ";
		}
		cout<<endl;
	}
	
	cout<<endl<<endl<<"Psi : "<<endl;
	for(int t=1; t<=T; t++){
		for(int j=1; j<=N; j++){
			cout<<Psi[t][j]<<" ";
		}
		cout<<endl;
	}
	*/

	/* termination */
	maxval = INT_MIN;
	maxarg = INT_MIN;
	for(int i=1; i<=N; i++) 
		if(Delta[T][i] > maxval) 
			maxval = Delta[T][i]  ,   maxarg = i;
	Pstar = maxval;
	//cout<<endl<<"Pstar : "<<Pstar<<endl;
	
	/******** State Seq. - Backtracking ********/
	Q[T] = maxarg;
	for(int t=T-1; t>=1; t--)
		Q[t] = Psi[t+1][Q[t+1]] ;

	/*
	cout<<endl<<endl<<"Obs. Sequence Taken -> O : "<<endl;
	for(int i=1; i<=T; i++){
		cout<<O[i]<<" ";
	}

	cout<<endl<<endl<<"State Sequence (Q) : "<<endl;
	for(int t=1; t<=T; t++){
		  cout<<Q[t]<<" ";
	}
	cout<<endl;
	*/
}


// Solution to Problem-3 : Xi, newPi, newA, newB
void prob3(){
	/******** Xi ********/
	//cout<<endl<<endl<<"Xi : "<<endl;
	for(int t=1; t<=T-1; t++){
		ld sum=0;
		for(int i=1; i<=N; i++){
			for(int j=1; j<=N; j++){ 
				Xi[t][i][j] = ( Alpha[t][i] * A[i][j] * B[j][O[t+1]] * Beta[t+1][j] ) ;
				sum += Xi[t][i][j];
			}
		}

		for(int i=1; i<=N; i++){
			for(int j=1; j<=N; j++){ 
				if(sum!=0) Xi[t][i][j] /= sum ;
				//cout<<Xi[t][i][j]<<" ";
			}
			//cout<<endl;
		}
		
		//cout<<endl;
	}


	/******** newPi ********/
	//cout<<endl<<"newPi : "<<endl;
	for(int i=1; i<=N; i++){
		newPi[i] = Gamma[1][i] ;
		//cout<<newPi[i]<<" ";
	}
	

	/******** newA ********/
	//cout<<endl<<endl<<"newA : "<<endl;
	ld num = 0, denom = 0;
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			num = 0, denom = 0;
			for(int t=1; t<=T-1; t++)
				num += Xi[t][i][j] ,   denom += Gamma[t][i]  ;
			if(denom!=0) newA[i][j] = num/denom ;
			else newA[i][j] = 0;
			//cout<<newA[i][j]<<" ";
		}
		// cout<<endl;
	}

	ld maxValue=-1, maxIndex=-1, sum=0;

	/******** newB ********/
	for(int j=1; j<=N; j++){
		maxValue=-1, maxIndex=-1, sum=0;
		for(int k=1; k<=M; k++){
			num = 0, denom = 0;
			for(int t=1; t<=T; t++){
				if(O[t] == k)
					num += Gamma[t][j] ; 
				denom += Gamma[t][j]  ;
			}
			if(denom==0){
				//cout<<"ERROR : DENOM = 0 !!"<<endl; 
				newB[j][k] = 0;
			}
			else newB[j][k] = num/denom ;

			if(newB[j][k] > maxValue) maxValue = newB[j][k], maxIndex = k;
			if(newB[j][k] < (ld)(pow((ld)10,(ld)-30))){
				sum += (ld)(pow((ld)10,(ld)-30)) - newB[j][k] ;
				newB[j][k] = (ld)(pow((ld)10,(ld)-30)) ;
			}
		}
		// cout<<maxValue<<" "<<maxIndex<<endl;
		newB[j][maxIndex] -= sum  ;
	}

	/*
	cout<<endl<<endl<<"newB : "<<endl;
	for(int j=1; j<=N; j++){
		for(int k=1; k<=M; k++){
			cout<<newB[j][k]<<" ";
		}
		cout<<endl;
	}
	*/

	/* updating vectors */
	A=newA;
	B=newB;
	Pi=newPi;
}

// creates the HMM Model for the current file
void createModel(){
	int itr = 1;         // iteration no.
	Pstar=0;
	prevPstar=-1;

	// while(prevPstar<Pstar){
	//while(prevPstar<Pstar && itr<=20){
	while(itr<=20){
		// cout<<endl<<endl<<endl<<"Iteration No : "<<iter-it<<endl;
		prob1();

		prob2();

		prob3();
		itr++;
	}
	cout<<"Iterations : "<<itr<<endl;
}

// dumps the current model of file in the appropriate folder under "Models" and further under the current iteration no.
void saveModel(string iterNo, string folderName, string prefixFile){
	string word, fileNameA, fileNameB, fileNamePi, dir;
	dir = "Models/ITER_" + iterNo + "/" + folderName;
	wstring wdir(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	dir = "Models/ITER_" + iterNo + "/" ;
	// cmd = "mkdir " + dir;
	// system(cmd.c_str());
	// cout<<"DIR : "<<dir<<endl;
	fileNameA =  dir + prefixFile + "_A.txt";
	fileNameB = dir + prefixFile + "_B.txt";
	fileNamePi = dir + prefixFile + "_Pi.txt";

	ofstream fout;
	fout.open(fileNameA);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fout<<A[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();

	fout.open(fileNameB);
        
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fout<<B[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();

	fout.open(fileNamePi);

	for(int i=1; i<=N; i++){
		fout<<Pi[i]<<" ";
	}
	fout<<endl;
	fout.close();
}

// loads the Observation Seq. to the "O" vector from the given "fileName"
void loadObsSeq(string fileName){
	//cout<<fileName<<endl;
	string word;
	ifstream fin(fileName);
	O.clear();
	O.push_back(-1);
	T=0;
	while(fin>>word){
		T++;
		O.push_back(stoi(word));
		// cout<<O.back()<<endl;
	}
	fin.close();

	cout<<"T : "<<T<<endl;

	Alpha = vector<vector<ld>> (T+1,vector<ld>(N+1));
	Beta = vector<vector<ld>> (T+1,vector<ld>(N+1,1)); // initialised to be 1
	Gamma = vector<vector<ld>> (T+1,vector<ld>(N+1));
	Delta = vector<vector<ld>> (T+1,vector<ld>(N+1));
	Psi = vector<vector<int>> (T+1,vector<int>(N+1,0)); // initialised to be 0
	Xi = vector<vector<vector<ld>>> (T+1,vector<vector<ld>>(N+1,vector<ld>(N+1)));
	Q = vector<int> (T+1);
}

// helps in generating the Avg. Model for A,B,Pi matrices for input to the next iteration
void genAvgModel(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			AvgA[i][j] += A[i][j];
		}
	}
       
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			AvgB[i][j] += B[i][j];
		}
	}

	for(int i=1; i<=N; i++){
		AvgPi[i] += Pi[i];
	}
}

// dumps the current avg model of the digit "k" in the Models folder
void updateFinalModel(string k){
	string fileNameA, fileNameB, fileNamePi;

	fileNameA = "Models/" + DEFAULT_ROLL + "_" + k + "_matrixA.txt";
	fileNameB = "Models/" + DEFAULT_ROLL + "_" + k + "_matrixB.txt";
	fileNamePi = "Models/" + DEFAULT_ROLL + "_" + k + "_matrixPi.txt";

	ofstream fout;
	fout.open(fileNameA);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fout<<AvgA[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
       
	fout.open(fileNameB);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fout<<AvgB[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();

	fout.open(fileNamePi);
	for(int i=1; i<=N; i++){
		fout<<AvgPi[i]<<" ";
	}
	fout<<endl;
	fout.close();
}

// dumps the current avg model of the digit "k" in the appropriate folder under "Models" and further under the current iteration no. 
void saveAvgModel(int numTrainfiles, string iterNo, string k){
	ld n = (ld)numTrainfiles;
	string fileNameA, fileNameB, fileNamePi;

	fileNameA = "Models/ITER_" + iterNo + "/" + ROLL + "_" + k + "_matrixA.txt";
	fileNameB = "Models/ITER_" + iterNo + "/" + ROLL + "_" + k + "_matrixB.txt";
	fileNamePi = "Models/ITER_" + iterNo + "/" + ROLL + "_" + k + "_matrixPi.txt";

	ofstream fout;
	fout.open(fileNameA);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			AvgA[i][j] /= n;
			fout<<AvgA[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
       
	fout.open(fileNameB);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			AvgB[i][j] /= n;
			fout<<AvgB[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();

	fout.open(fileNamePi);
	for(int i=1; i<=N; i++){
		AvgPi[i] /= n;
		fout<<AvgPi[i]<<" ";
	}
	fout<<endl;
	fout.close();
}

// selects the best model for a digit
void selectBestModel(){
	if(Pstar > BestPstar){
		BestPstar = Pstar;
		BestA = A;
		BestB = B;
		BestPi = Pi;
	}
}

// dumps the best model for digit "k" in the "Models" folder
void saveBestModel(string k){
	string fileNameA, fileNameB, fileNamePi;

	fileNameA = "Models/" + ROLL + "_" + k + "_BestMatrixA.txt";
	fileNameB = "Models/" + ROLL + "_" + k + "_BestMatrixB.txt";
	fileNamePi = "Models/" + ROLL + "_" + k + "_BestMatrixPi.txt";

	ofstream fout;
	fout.open(fileNameA);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fout<<BestA[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
       
	fout.open(fileNameB);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fout<<BestB[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();

	fout.open(fileNamePi);
	for(int i=1; i<=N; i++){
		fout<<BestPi[i]<<" ";
	}
	fout<<endl;
	fout.close();
}

// generates an avg model for each digit w.r.t. all the training files and dumps it in the appropriate folder under "Models" directory
void GenerateModels(int numTrainfiles){
	int ITER = 1;        // model training iteration no.
	string fileName, folderName, word, prefixFile, dir;

	dir = "Models";
	wstring wdir(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	string newname = CUR_WORD;
	vector<string> newfile;
	for(int i=1; i<=numTrainfiles; i++){
		stringstream ss;
		ss << i;
		newfile.push_back(ss.str());
	}

	bool flg=0;

	while(ITER <= TOTAL_ITER){
		stringstream ss;
		ss << ITER;

		dir = "Models/ITER_" + ss.str();
		wdir.clear();
		wdir.assign(dir.begin(), dir.end());
		if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
		else cout<<"False "<<endl;

		if(ITER == TOTAL_ITER) flg=1;
		for(int k=0; k < 1; k++){
			folderName = ROLL + "_" + newname ;
			AvgPi = vector<ld> (N+1,0); 
			AvgA = vector<vector<ld>> (N+1,vector<ld>(N+1,0));
			AvgB = vector<vector<ld>> (N+1,vector<ld>(M+1,0));
			
			if(flg){
				BestPi = vector<ld> (N+1); 
				BestA = vector<vector<ld>> (N+1,vector<ld>(N+1));
				BestB = vector<vector<ld>> (N+1,vector<ld>(M+1));
				BestPstar = -1;
			}

			for(int fid=0; fid<numTrainfiles; fid++){
				cout<<endl<<"ITER : "<<ITER<<",  k : "<<k<<",  file_no. : "<<fid+1<<endl;
				initialiseModel(ss.str(), newname);  // populating the initial model A, B, Pi matrices for the digit 'k'
				prefixFile =  folderName + "/" + ROLL + "_" + newname + "_" + newfile[fid];
				fileName = prefixFile + ".txt";
				cout<<prefixFile<<endl;
				loadObsSeq("ObsvSeq/" + fileName);
				createModel();
				saveModel(ss.str(), folderName, prefixFile);

				if(flg) selectBestModel();
				genAvgModel();
			}

			if(flg) saveBestModel(newname);
			saveAvgModel(numTrainfiles, ss.str(), newname);
			updateFinalModel(newname);
		} 

		ITER++;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////// TestModels ///////////////////////////////////////////////////////////////
// predicted, actual  : stores predicted & actual values resp. per digit
vector<string> predicted, actual;

// loads the final Model to be used for Testing for the digit "k"
void loadFinalModel(string k){
	string fileNameA, fileNameB, fileNamePi, word;
	/*
	fileNameA = "Models/" + ROLL + "_" + k + "_BestMatrixA.txt";
	fileNameB = "Models/" + ROLL + "_" + k + "_BestMatrixB.txt";
	fileNamePi = "Models/" + ROLL + "_" + k + "_BestMatrixPi.txt";
	*/

	fileNameA = "Models/" + ROLL + "_" + k + "_matrixA.txt";
	fileNameB = "Models/" + ROLL + "_" + k + "_matrixB.txt";
	fileNamePi = "Models/" + ROLL + "_" + k + "_matrixPi.txt";
	cout<<fileNameA<<endl;
	ifstream fin;
	fin.open(fileNameA);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fin >> word;
			A[i][j] = stold(word);
		}
	}
	fin.close();

	fin.open(fileNameB);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fin >> word;
			B[i][j] = stold(word);
		}
	}
	fin.close();

	fin.open(fileNamePi);
	for(int i=1; i<=N; i++){
		fin >> word;
		Pi[i] = stold(word);
	}
	fin.close();
}

// predicts the digit spoken in the current file by taking that digit's model having maximum P(O/lambda) w.r.t. the observation seq. of the current file 
// and stores the response( predicteddigit) in the "predicted" vector
void PredictDigit(string fileName, string word){
	
	ld maxP = -1;
	int maxPIndx = -1;
	loadObsSeq(fileName);
	vector<string> wordArr;
	if(word == "BEGIN"){
		wordArr.push_back("BEGIN");
	}
	else if(word == "POSITION"){
		wordArr.push_back("PAUSE");
		for(int i=0; i<9; i++){
			stringstream ss;
			ss << i;
			wordArr.push_back(ss.str());
		}
	}
	else if(word == "CONFIRM"){
		wordArr.push_back("YES");
		wordArr.push_back("NO");
		//wordArr.push_back("PAUSE");
	}
	else if(word == "PAUSE_MENU"){
		wordArr.push_back("RESUME");
		//wordArr.push_back("RESTART");
		wordArr.push_back("BEGIN");
		wordArr.push_back("QUIT");
	}
	else{
		for(int i=0; i<10; i++){
			stringstream ss;
			ss << i;
			wordArr.push_back(ss.str());
		}
	}

	for(int i=0;i<wordArr.size();i++){
		loadFinalModel(wordArr[i]);
		Alpha = vector<vector<ld>> (T+1,vector<ld>(N+1));
		Beta = vector<vector<ld>> (T+1,vector<ld>(N+1,1)); // initialised to be 1
		prob1();
		if(ProbOgivenLambda > maxP) maxP=ProbOgivenLambda, maxPIndx=i;
	}
	//cout<<"maxP : "<<maxP<<endl;
	predicted.push_back(wordArr[maxPIndx]);

	return;
}

// processes the testing files and stores the predicted spoken digit for each testing file per digit
void TestCompute(int numTestfiles, int numTrainfiles){
	string digit = CUR_WORD;
	string fileName;
	for(int j=0;j<1;j++){
		for(int k=numTrainfiles+1;k<=numTrainfiles+numTestfiles;k++){
			stringstream ss;
			ss << k;
			cout<<endl<<"Digit being tested : "<<digit<<"  file_num : "<<k<<endl;
			fileName = "ObsvSeq/" + ROLL + "_" + digit + "/" + ROLL + "_" + digit + "_" + ss.str() + ".txt" ;
			PredictDigit(fileName, "");
		}
	}
}

// test the training files for the correctness of the predicted digits w.r.t actual spoken ones and report the accuracy
void TestModels(int numTestfiles, int numTrainfiles){
	int numCorrect = 0;
	for(int i=0;i<10;i++){
		stringstream ss;
		ss << i;
		for(int j=0;j<numTestfiles;j++){
			actual.push_back(ss.str());
		}
	}

	TestCompute(numTestfiles, numTrainfiles);

	cout<<endl<<"Actual"<<endl;
	for(int i=0;i<actual.size();i++){
		cout<<actual[i]<<" ";
	}
	cout<<endl;

	cout<<endl<<"Predicted"<<endl;
	for(int i=0;i<predicted.size();i++){
		cout<<predicted[i]<<" ";
	}
	cout<<endl;

	// comparing actual and predicted digits
	for(int i=0;i<actual.size();i++){
		if(actual[i] == predicted[i]) numCorrect++;
	}
	if(actual.size() != 0) cout<<endl<<"Accuracy : "<< (100.0 * (ld)(numCorrect)) /(ld)(actual.size())<<endl;
	else cout<<"Test files : 0"<<endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////// TestRecCompute ///////////////////////////////////////////////////////////
// processes the testing file and stores the predicted spoken digit for the current recording of the digit
void TestRecCompute(string codeBookFileName, string word){

	string fileName, new_func, new_name, time;

	// recording_module 
	time = "3 ";  // no. of seconds for recording
	new_name = "0_rec_1"; // new name of files (both .wav and .txt)
	//new_func = "Recording_Module.exe " + time + new_name + ".wav " + new_name + ".txt"; 
	new_func = "Recording_Module.exe " + time + new_name + ".wav " + new_name + ".txt & exit"; 
	system(new_func.c_str());
	//WinExec(new_func.c_str(),0);

	fileName = new_name + ".txt";

	cout<<"REC. FILE NAME : "<<fileName<<endl;

	R_file.clear();
	A_file.clear(); 
	C_file.clear();
	processData(fileName, fileName);
	GenerateObsvSeq(0, codeBookFileName, 1, fileName);
	PredictDigit("ObsvSeq/"+fileName, word);
	 
	return;		
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<char> matrix(9, '-');
char player[2] = {'O','X'};
int flg=1;
int hours = 0; 
int minutes = 0; 
int seconds = 0; 
  
// function to display the timer 
void displayMatrix() 
{ 
    // system call to clear the screen 
    //system("clear"); 
    //system("CLS"); 
    cout << setfill(' ') << setw(55) << "      MATRIX GRID      \n"; 

    cout << setfill(' ') << setw(46) << " -----------------\n";

    cout << setfill(' ') << setw(29); 
    cout << "| " << setw(2) << matrix[0] << "  | "; 
    cout  << setw(2) << matrix[1] << "  | "; 
    cout  << setw(2) << matrix[2] << "  |" << endl; 
    cout << setfill(' ') << setw(46) << " -----------------\n";

    cout << setfill(' ') << setw(29); 
    cout << "| " << setw(2) << matrix[3] << "  | "; 
    cout  << setw(2) << matrix[4] << "  | "; 
    cout << setw(2) << matrix[5] << "  |" << endl; 
    cout << setfill(' ') << setw(46) << " -----------------\n";

    cout << setfill(' ') << setw(29); 
    cout << "| "  << setw(2) << matrix[6] << "  | "; 
    cout  << setw(2) << matrix[7] << "  | "; 
    cout  << setw(2) << matrix[8] << "  |" << endl; 
    cout << setfill(' ') << setw(46) << " -----------------\n";
} 

void displayClock() 
{ 
    // system call to clear the screen 
    //system("clear"); 
    system("CLS"); 
    cout << setfill(' ') << setw(62) << "         TIME LEFT TILL SPEECH START       \n"; 
    cout << setfill(' ') << setw(55) << " --------------------------\n"; 
    cout << setfill(' ') << setw(29); 
    cout << "| " << setfill('0') << setw(2) << hours << " hrs | "; 
    cout << setfill('0') << setw(2) << minutes << " min | "; 
    cout << setfill('0') << setw(2) << seconds << " sec |" << endl; 
    cout << setfill(' ') << setw(55) << " --------------------------\n"; 

}

void timer() 
{ 
    // infinte loop because timer will keep  
    // counting. To kill the process press 
    // Ctrl+D. If it does not work ask 
    // ubuntu for other ways. 
    while (true) { 
          
        // display the timer 
        displayClock(); 
  
        // sleep system call to sleep  
        // for 1 second 
        Sleep(1000); 
  
        // increment seconds 
        seconds++; 
  
        // if seconds reaches 60 
        if (seconds == 60) { 
          
            // increment minutes 
            minutes++; 
  
            // if minutes reaches 60 
            if (minutes == 60) { 
          
                // increment hours 
                hours++; 
                minutes = 0; 
            } 
            seconds = 0; 
        } 
    } 
} 

///////////////////////////////////////////////////////////////////////////////////////////////////////


void detectBegin(){
	while(1){
		for(int i=3; i>0; i--){
			seconds=i;
			Sleep(1000);
			displayClock();
			cout<<"At prompt, PLEASE speak \"BEGIN\" "<<endl;
		}
		seconds=0;
		displayClock();
		cout<<"At prompt, PLEASE speak \"BEGIN\" "<<endl;
		Sleep(2000);

		TestRecCompute(genCodeBookfile,"BEGIN");
		string word = predicted[0];
		cout<<word<<endl;
		cout<<"Word Predicted : "<<word<<endl<<endl;
		predicted.clear();
		Sleep(2000);
		if(word=="BEGIN") break;
	}
}

string detectPosition(){
	string word;
	while(1){
		int last=0;
		for(int i=3; i>0; i--){
			seconds=i;
			Sleep(1000);
			displayClock();
			displayMatrix();
			cout<<"At prompt, PLAYER-"<<player[flg]<<"  ROLL: "<<ROLL<<" :: PLEASE SPEAK one of the DIGITS as positions of the matrix :-  {";
			for(int j=0;j<9;j++)
				if(matrix[j]=='-')
					last=j;
			for(int j=0;j<9;j++)
				if((matrix[j]=='-')&&(j!=last))
					cout<<j<<",";
			cout<<last<<"} \n";
		}
		seconds=0;
		displayClock();
		displayMatrix();
		cout<<"At prompt, PLAYER-"<<player[flg]<<"  ROLL: "<<ROLL<<" :: PLEASE SPEAK one of the DIGITS as positions of the matrix :-  {";
		for(int j=0;j<9;j++)
			if(matrix[j]=='-')
				last=j;
		for(int j=0;j<9;j++)
			if((matrix[j]=='-')&&(j!=last))
				cout<<j<<",";
		cout<<last<<"} \n";
		Sleep(2000);

		TestRecCompute(genCodeBookfile,"POSITION");
		word = predicted[0];
		cout<<word<<endl;
		cout<<"Word Predicted : "<<word<<endl<<endl;
		predicted.clear();
		Sleep(2000);
		if(word=="PAUSE" || word=="0" || word=="1" || word=="2" || word=="3" || word=="4" || word=="5" || word=="6" || word=="7" || word=="8" || word=="9") break;
	}
	return word;
}

string detectConfirmation(){
	string word;
	while(1){
		for(int i=3; i>0; i--){
			seconds=i;
			Sleep(1000);
			displayClock();
			displayMatrix();
			cout<<"At prompt, PLAYER-"<<player[flg]<<"  ROLL: "<<ROLL<<" :: PLEASE CONFIRM whether matrix updation is correct or not by speaking either :- \"YES\" OR \"NO\" "<<endl;
		}

		seconds=0;
		displayClock();
		displayMatrix();
		cout<<"At prompt, PLAYER-"<<player[flg]<<"  ROLL: "<<ROLL<<" :: PLEASE CONFIRM whether matrix updation is correct or not by speaking either :- \"YES\" OR \"NO\" "<<endl;
		Sleep(2000);

		TestRecCompute(genCodeBookfile,"CONFIRM");
		word = predicted[0];
		cout<<word<<endl;
		cout<<"Word Predicted : "<<word<<endl<<endl;
		predicted.clear();
		Sleep(2000);
		if(word=="PAUSE" || word=="YES" || word=="NO") break;
	}
	return word;
}

string detectOption(){
	string word;
	while(1){
		for(int i=3; i>0; i--){
			seconds=i;
			Sleep(1000);
			displayClock();
			//displayMatrix();
			cout<<"At prompt, PLAYER-"<<player[flg]<<"  ROLL: "<<ROLL<<" :: PLEASE SELECT one of the following options by speaking them :- \"RESUME\" OR \"BEGIN(RESTART)\" OR \"QUIT\" "<<endl;
		}

		seconds=0;
		displayClock();
		displayMatrix();
		cout<<"At prompt, PLAYER-"<<player[flg]<<"  ROLL: "<<ROLL<<" :: PLEASE SELECT one of the following options by speaking them :- \"RESUME\" OR \"BEGIN(RESTART)\" OR \"QUIT\" "<<endl;
		Sleep(2000);

		TestRecCompute(genCodeBookfile,"PAUSE_MENU");
		word = predicted[0];
		cout<<word<<endl;
		cout<<"Word Predicted : "<<word<<endl<<endl;
		predicted.clear();
		Sleep(2000);

		if(word=="RESUME" || word=="BEGIN" || word=="QUIT") break;
		//if(word=="RESUME" || word=="RESTART" || word=="QUIT") break;
	}
	return word;
}

int solve(){
	int mat[9];
	int flag=0;
	for(int i=0;i<9;i++){
			if(matrix[i] == 'X')
				mat[i]=1;
			else if(matrix[i] == 'O')
				mat[i]=-1;
			else{
				mat[i]=0;
				flag=1;
			}
	}
	int sum1,sum2;
	for(int i=0;i<3;i++){
		sum1=0;
		sum2=0;
		for(int j=0;j<3;j++){
			sum1+=mat[3*i+j];
			sum2+=mat[i+(3*j)];
		}
		if((sum1==3)||(sum2==3))
			return 1;
		if((sum1==-3)||(sum2==-3))
			return 2;
	}
	sum1=mat[0]+mat[4]+mat[8];
	sum2=mat[2]+mat[4]+mat[6];
	if((sum1==3)||(sum2==3))
		return 1;
	if((sum1==-3)||(sum2==-3))
		return 2;
	if(!flag)
		return -1;
	return 0;
}

void GAME() 
{ 
    // infinte loop because timer will keep  
    // counting. To kill the process press 
    // Ctrl+D. If it does not work ask 
    // ubuntu for other ways. 
   Sleep(2000);
   ROLL = ORIGINAL_ROLL;
   detectBegin();
   Sleep(2000);
   matrix = vector<char> (9, '-');
   displayMatrix();

   Sleep(2000);

   flg=1;  // flg: 1 --> X,   flg: 0 ---> O,    player-1 (X) starts
   int flag;
   
   string detectedWord, detectConf, detectedOption;
   while(1){
	   flag=0;
	   if(flg) ROLL = rollPlayer1;
	   else ROLL = rollPlayer2;

	   detectedWord = detectPosition();
	   Sleep(2000);
	   if(detectedWord == "PAUSE"){
		   detectedOption = detectOption();
		   Sleep(2000);
		   if(detectedOption == "RESUME"){
			   cout<<"DETECTED::RESUME :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : PLEASE SPEAK THE POSITION AGAIN"<<endl;
			   Sleep(2000);
			   continue;
		   }
		   //else if(detectedOption == "RESTART"){
		   else if(detectedOption == "BEGIN"){
			   //cout<<"DETECTED::RESTART :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : WANTS TO RESTART AGAIN"<<endl;
			   cout<<"DETECTED::BEGIN :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : WANTS TO RESTART AGAIN"<<endl;
			   cout<<"RESTARTING..."<<endl;
			   GAME();
			   Sleep(2000);
			   break;
		   }
		   else if(detectedOption == "QUIT"){
			   cout<<"DETECTED::QUIT :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : WANTS TO QUIT THE GAME"<<endl;
			   cout<<"QUITTING..."<<endl;
			   Sleep(2000);
			   break;
		   }
	   }
	   else{
		   //ask for confirmation
		   detectConf = detectConfirmation();
		   Sleep(2000);
		   if(detectConf == "YES"){
			   if(detectedWord[0]<'0' || detectedWord[0]>'8'){
				   cout<<"INVALID CHOICE, PLEASE TRY AGAIN\n";
				   flag=1;
			   }
			   else if((matrix[stoi(detectedWord)]=='O')||(matrix[stoi(detectedWord)]=='X')){
				   cout<<"WRONG CHOICE, PLEASE TRY AGAIN\n";
				   flag=1;
			   }
			   else{
				   matrix[stoi(detectedWord)] = player[flg];
				   // check for repeated position and winning/draw condition
				   // if win/draw, same as QUIT, with a prompt of which player wins finally
				   Sleep(2000);
				   int x = solve();
				   //int x=0;
				   cout<<endl;
				   if(x<0){
					   cout<<"MATCH DRAWN!"<<endl;
					   Sleep(5000);
					   break;
				   }
				   else if(x>0){
					   cout<<"PLAYER "<<x<<" WON!"<<endl;
					   Sleep(5000);
					   break;
				   }
			   }
		   }
		   else if(detectConf == "NO"){
			   cout<<"DETECTED::NO :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : PLEASE SPEAK THE POSITION AGAIN"<<endl;
			   Sleep(2000);
			   continue;
		   }
		   else if(detectConf == "PAUSE"){

			   cout<<"DETECTED::PAUSE"<<endl;
			   detectedOption = detectOption();
			   Sleep(2000);
			   if(detectedOption == "RESUME"){
				   cout<<"DETECTED::RESUME :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : PLEASE SPEAK THE POSITION AGAIN"<<endl;
				   Sleep(2000);
				   continue;
			   }
			   //else if(detectedOption == "RESTART"){
			   else if(detectedOption == "BEGIN"){
				   //cout<<"DETECTED::RESTART :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : WANTS TO RESTART AGAIN"<<endl;
				   cout<<"DETECTED::BEGIN :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : WANTS TO RESTART AGAIN"<<endl;
				   cout<<"RESTARTING..."<<endl;
				   GAME();
				   Sleep(2000);
				   break;
			   }
			   else if(detectedOption == "QUIT"){
				   cout<<"DETECTED::QUIT :: PLAYER "<<flg<<" ("<<player[flg]<<")  ROLL: "<<ROLL<<"  : WANTS TO QUIT THE GAME"<<endl;
				   cout<<"QUITTING..."<<endl;
				   Sleep(2000);
				   break;
			   }
		   }
	   }
	   if(!flag)
			flg = 1 - flg;
   }
   
} 

///////////////////////////////////////////////////////
void trainableModule(string roll){
	string dir;
	ROLL = roll;
	Recordings_Dir = ROLL+"_Recordings/";

	dir = Recordings_Dir;
	wstring wdir(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	//*
	//system("Recording_Module.exe 3 test.wav test.txt");
	//char digit[10]={'0','1','2','3','4','5','6','7','8','9'};
	string foldername,filename;
	for(int i=0;i<10;i++){
		cout<<"SPEAK DIGIT: "<<i<<"  :: "<<TRAINABLE_MODULE_REC<<" times\n";
		foldername = Recordings_Dir+roll+"_"+to_string((long long)i);
		dir = foldername;
		wstring wdir2(dir.begin(), dir.end());
		if (CreateDirectory(wdir2.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
		else cout<<"False "<<endl;
		//string createfolder="mkdir "+foldername;
		//system(createfolder.c_str());
		for(int j=0;j<TRAINABLE_MODULE_REC;j++){
			cout<<"SPEAK DIGIT: "<<i<<"  :: "<<TRAINABLE_MODULE_REC-j<<" more times\n";
			filename = roll+"_"+to_string((long long)i)+"_"+to_string((long long)j);
			//cout<<filename<<"\n";
			string f1 = filename+".wav";
			string f2 = filename+".txt";
			cout<<f1<<" "<<f2<<"\n";
			string command = "Recording_Module.exe 1 " + foldername+"\\"+f1 + " " + foldername+"\\"+f2;
			system(command.c_str());
		}
	}
	string s[7]={"BEGIN","NO","PAUSE","QUIT","RESTART","RESUME","YES"};
	for(int i=0;i<sizeof(s)/sizeof(s[0]);i++){
		cout<<"SPEAK WORD: "<<s[i]<<"  :: "<<TRAINABLE_MODULE_REC<<" times\n";
		foldername = Recordings_Dir+roll+"_"+s[i];
		dir = foldername;
		wstring wdir2(dir.begin(), dir.end());
		if (CreateDirectory(wdir2.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
		else cout<<"False "<<endl;

		//string createfolder="mkdir "+foldername;
		//system(createfolder.c_str());
		for(int j=0;j<TRAINABLE_MODULE_REC;j++){
			cout<<"SPEAK WORD: "<<s[i]<<"  :: "<<TRAINABLE_MODULE_REC-j<<" more times\n";
			string filename=roll+"_"+s[i]+"_"+to_string((long long)j);
			string f1 = filename+".wav";
			string f2 = filename+".txt";
			cout<<f1<<" "<<f2<<"\n";
			string command = "Recording_Module.exe 1 " + foldername+"\\"+f1 + " " + foldername+"\\"+f2;
			system(command.c_str());
		}
	}
	//*/
}

void createNewModel(string roll){
	trainableModule(roll);
}

int selectPlayerOption(int playernum){
	int opt=-1;
	while(1){
		cout<<"WHICH PLAYER DO YOU WANT PLAYER \""<<playernum<<"\" TO PLAY AS?"<<endl;
		cout<<"Select one of the option numbers by entering it..."<<endl<<endl;
		cout<<1<<" : PLAYER A"<<endl;
		cout<<2<<" : PLAYER B"<<endl;
		cout<<3<<" : CUSTOM PLAYER"<<endl;
		cout<<4<<" : QUIT"<<endl;

		cin>>opt;
		if(opt<1 || opt>4) cout<<"PLEASE ENTER VALID OPTION NO. AGAIN!"<<endl<<endl;
		else break;
	}
	
	return opt;
}

int selectMainOption(){
	int opt=-1;
	while(1){
		cout<<"WHICH OPTION DO YOU WANT ?"<<endl;
		cout<<"Select one of the option numbers by entering it..."<<endl<<endl;
		cout<<1<<" : PLAY GAME"<<endl;
		cout<<2<<" : TRAIN A WORD"<<endl;
		cout<<3<<" : QUIT"<<endl;

		cin>>opt;
		if(opt<1 || opt>3) cout<<"PLEASE ENTER VALID OPTION NO. AGAIN!"<<endl<<endl;
		else break;
	}
	
	return opt;
}

int selectOption(){
	int opt=-1;
	while(1){
		cout<<"WHICH DIGIT/WORD DO YOU WANT TO TRAIN ?"<<endl;
		cout<<"Select one of the option numbers by entering it..."<<endl<<endl;
		for(int i=0; i<10; i++){
			cout<<i<<" : DIGIT \""<<i<<"\" "<<endl;
		}
		string words[] = {"BEGIN","NO","PAUSE","QUIT","RESTART","RESUME","YES"};
		for(int i=10; i<17; i++){
			cout<<i<<" : WORD \""<<words[i-10]<<"\" "<<endl;
		}
		cout<<17<<" : QUIT"<<endl;

		cin>>opt;
		if(opt<0 || opt>17) cout<<"PLEASE ENTER VALID OPTION NO. AGAIN!"<<endl<<endl;
		else break;
	}
	
	return opt;
}

void trainWord(int opt){
	string roll = DEFAULT_ROLL;

	string word = WORDS[opt];
	CUR_WORD = word;

	string dir;
	ROLL = roll;
	Recordings_Dir = ROLL+"_Recordings/";

	dir = Recordings_Dir;
	wstring wdir(dir.begin(), dir.end());
	if (CreateDirectory(wdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	//*
	//system("Recording_Module.exe 3 test.wav test.txt");
	//char digit[10]={'0','1','2','3','4','5','6','7','8','9'};
	string foldername,filename;

	cout<<"SPEAK WORD: "<<word<<"  :: "<<TRAINABLE_MODULE_REC<<" times\n";
	foldername = Recordings_Dir+roll+"_"+word;
	dir = foldername;
	wstring wdir2(dir.begin(), dir.end());
	if (CreateDirectory(wdir2.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()); //cout<<"TRUE "<<endl;
	else cout<<"False "<<endl;

	numTrainfiles=TRAINABLE_MODULE_REC;
	numTestfiles=0;

	
	for(int j=1;j<=TRAINABLE_MODULE_REC;j++){
		cout<<"SPEAK WORD: "<<word<<"  :: "<<TRAINABLE_MODULE_REC+1-j<<" more times\n";
		string filename=roll+"_"+word+"_"+to_string((long long)j);
		string f1 = filename+".wav";
		string f2 = filename+".txt";
		cout<<f1<<" "<<f2<<"\n";
		string command = "Recording_Module.exe 3 " + foldername+"\\"+f1 + " " + foldername+"\\"+f2;
		system(command.c_str());
	}
	
	Compute_Ci(numTrainfiles+numTestfiles);
	
	// generates the CodeBook using the Universe Ci vectors in the given file with name as : "genUniversefile" 
	//genCodeBook(genUniversefile);

	// Generates the Obs. seq for all the recorded files and dumps each one of them in the appropriate folder under "ObsvSeq" folder
	GenerateObsvSeq(numTrainfiles+numTestfiles, genCodeBookfile, 0, "");
	

	// generates an avg model for each digit w.r.t. all the training files and dumps it in the appropriate folder under "Models" directory
	GenerateModels(numTrainfiles);

}

///////////////////////////////////////////////////////
int main()
{
	// rec_flg = 0  --> no recording, accuracy of test data is computed
	// rec_flg = 1  --> recording on, accuracy of the spoken digit is to be computed
	bool rec_flg = 1; 

	cout << setprecision(30);
	//numTrainfiles = min(numTrainfiles, 29);
	//numTestfiles = min(numTestfiles, 30-numTrainfiles);
	//cout<<"numTrainfiles : "<<numTrainfiles<<" ;   numTestfiles : "<<numTestfiles<<endl<<endl;
	
	//string codeBookFileName = "original_codebook.txt";  // select the given Codebook file
	string codeBookFileName = genCodeBookfile; // select the LBG-generated Codebook file
	/*
	if(rec_flg==0){
		// processes the training & testing files, saves Ci vectors (& Trimmed Recordings) for each file in appropriate folder under "CiValues" folder (& "TrimmedRecordings" folder resp.), 
		// and also dumps the universe Ci vector 
		
		Compute_Ci(numTrainfiles+numTestfiles);
	
		// generates the CodeBook using the Universe Ci vectors in the given file with name as : "genUniversefile" 
		genCodeBook(genUniversefile);

		// Generates the Obs. seq for all the recorded files and dumps each one of them in the appropriate folder under "ObsvSeq" folder
		GenerateObsvSeq(numTrainfiles+numTestfiles, codeBookFileName, 0, "");
		
		// generates an avg model for each digit w.r.t. all the training files and dumps it in the appropriate folder under "Models" directory
		GenerateModels(numTrainfiles);
		
		// test the training files for the correctness of the predicted digits w.r.t actual spoken ones and report the accuracy
		TestModels(numTestfiles,numTrainfiles);
	}
	else{
		TestRecCompute(codeBookFileName);
		cout<<predicted[0]<<endl;
		cout<<"Digit Predicted : "<<predicted[0]<<endl<<endl;
	}
	*/
	//timer();

	string ans, pathname;
	int opt;
	opt = selectMainOption();
	if(opt == 3){
		cout<<"OPTION SELECTED::QUIT"<<endl<<"QUITTING..."<<endl;
		return 0;
	}
	else if(opt == 2){
		opt = selectOption();
	
		if(opt == 17){
			cout<<"OPTION SELECTED::QUIT"<<endl<<"QUITTING..."<<endl;
			return 0;
		}
		else{
			trainWord(opt);
			system("pause");
			system("CLS");
			opt = main();
			return 0;
		}
	}
	
	opt = selectPlayerOption(1);

	if(opt == 4){
		cout<<"OPTION SELECTED::QUIT"<<endl<<"QUITTING..."<<endl;
		return 0;
	}
	else{
		if(opt==1) rollPlayer1 = ROLL1;
		else if(opt==2) rollPlayer1 = ROLL2;
		else rollPlayer1 = DEFAULT_ROLL;

		opt = selectPlayerOption(2);
		if(opt == 4){
			cout<<"OPTION SELECTED::QUIT"<<endl<<"QUITTING..."<<endl;
			return 0;
		}
		else{
			if(opt==1) rollPlayer2 = ROLL1;
			else if(opt==2) rollPlayer2 = ROLL2;
			else rollPlayer2 = DEFAULT_ROLL;

			system("CLS");
		}
	}

	/*
	pathname = rollPlayer1+"_Recordings";

	if( stat( pathname.c_str(), &info ) != 0 ){
		cout<<"Cannot access "<<pathname<<"\n";
		cout<<"Please Record the words and digits first!!"<<endl;
		cout<<"Do you want to record now ? (Y/N)"<<endl;
		cin>>ans;
		transform(ans.begin(), ans.end(), ans.begin(), ::tolower);

		if(ans=="y" || ans=="yes"){
			createNewModel(rollPlayer1);
		}
		else return 0;
	}
	else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
		cout<<"Good! "<<pathname<<" directory exists."<<"\n";
	else{
		cout<<pathname<<" is NOT a directory."<<"\n";
		cout<<"Please Record the words and digits first!!"<<endl;
		cout<<"Do you want to record now ? (Y/N)"<<endl;
		cin>>ans;
		transform(ans.begin(), ans.end(), ans.begin(), ::tolower);

		if(ans=="y" || ans=="yes"){
			createNewModel(rollPlayer1);
		}
		else return 0;
	}


	cout<<endl<<endl<<"PLEASE ENTER 2nd PLAYER'S ROLL NUMBER:"<<"\n";
	cin>>rollPlayer1;
	pathname = rollPlayer1+"_Recordings";

	if( stat( pathname.c_str(), &info ) != 0 ){
		cout<<"Cannot access "<<pathname<<"\n";
		cout<<"Please Record the words and digits first!!"<<endl;
		cout<<"Do you want to record now ? (Y/N)"<<endl;
		cin>>ans;
		transform(ans.begin(), ans.end(), ans.begin(), ::tolower);

		if(ans=="y" || ans=="yes"){
			createNewModel(rollPlayer1);
		}
		else return 0;
	}
	else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
		cout<<"Good! "<<pathname<<" directory exists."<<"\n";
	else{
		cout<<pathname<<" is NOT a directory."<<"\n";
		cout<<"Please Record the words and digits first!!"<<endl;
		cout<<"Do you want to record now ? (Y/N)"<<endl;
		cin>>ans;
		transform(ans.begin(), ans.end(), ans.begin(), ::tolower);

		if(ans=="y" || ans=="yes"){
			createNewModel(rollPlayer1);
		}
		else return 0;
	}
	*/
	GAME();
	




	return 0;
}

