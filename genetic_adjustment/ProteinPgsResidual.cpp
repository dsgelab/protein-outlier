#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <armadillo>


#define NAcode -99999
using namespace std;
char buffer[100000];

typedef struct Sample {
	vector<double> PGS;
	vector<double> Pheno;
	vector<double> Covar;
	vector<double> StandRes;
	double AggStandRes = NAcode;
} IndInfo;
std::unordered_map<string, IndInfo> SampleList;

typedef struct Trait {
	int Idx;
	int HasPheno = 0;
	double R2 = 0.0;
	double adjR2 = 0.0;
	double Corr = 0.0;
} TraitInfo;
std::unordered_map<string, TraitInfo> TraitList;
int UseCovar;
int nTrait;
int nCovar;
long int nSample;
double R2thres;

string InGenoFile;
string InPhenoFile;
string InCovarFile;
string OutFile;

arma::mat GeneMat;
arma::mat PhenoMat;
arma::mat CovarMat;


void ReadParam(const char *ParFile) {
	FILE *Par;
	char *tok; char *p;
	R2thres = 0;
	UseCovar = 1;
	Par = fopen(ParFile,"r");
	// Set default values
	
	if (Par == NULL) {
	    printf("Cannot open parameter file.\n");
	    exit(0);
	}
	else {
		while (fgets(buffer, sizeof(buffer), Par) != NULL) {
			p = buffer;
			tok = strtok_r(p, " \t", &p);
	    	if (tok != NULL) {
	    		if (strcmp(tok, "GenoFileName") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			InGenoFile = tok;
				}
				else if (strcmp(tok, "PhenoFileName") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			InPhenoFile = tok;
				}
	    		else if (strcmp(tok, "CovarFileName") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			InCovarFile = tok;
				}
				else if (strcmp(tok, "R2cutoff") == 0) {
					tok = strtok_r(p, " \t\n", &p);
					R2thres = atof(tok);
					if (R2thres < 0) {
	    				printf("R2 cutoff for phenotype inclusion should be non-negative.\n");
						exit(0);
	    			}

				}
	    		else if (strcmp(tok, "OutFileName") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			OutFile = tok;
	    		}
	    	}
		}
	}

	if (InGenoFile.empty() ||  InPhenoFile.empty()) {
		printf("Missing genetic risk file or phenotype file.\n");
		exit(0);
	}
	if (InCovarFile.empty()) {
		printf("No covariate file input. Will not include covariate in regressions.\n");
		UseCovar = 0;
		nCovar = 0;
	}
	if (OutFile.empty()) {
		OutFile = "Output";
	}
	printf("Read parameters, done.\n");
	printf("Results will write to file %s.\n", OutFile.c_str());
	fclose(Par);
}



void ReadGenoFile(const char *FileName) {
	char *tok; char *p;
	long int i = 0;
	nTrait = 0;
	IndInfo tmpInd;
	string tmpIndIdx;
	double tmpVal;

	FILE *GenoFile;
	GenoFile = fopen(FileName, "r");

	if (GenoFile == NULL) {
		printf("Cannot open genetic risk file %s.\n", FileName);
		exit(0);
	}
	else {
		i = 0;
		fgets(buffer, sizeof(buffer), GenoFile); // header line
		p = buffer;
		tok = strtok_r(p, " ,\t\n", &p);
		while ((tok = strtok_r(p, " ,\t\n", &p))) {
			TraitList[tok].Idx = nTrait;
			nTrait++;
		}
		vector<double> tmpStandRes(nTrait, NAcode);
		
		i = 0; // Individual counter
		while (fgets(buffer, sizeof(buffer), GenoFile) != NULL) {
			p = buffer;
			tok = strtok_r(p, " ,\t\n", &p); //ID
			tmpIndIdx = tok;
			while ((tok = strtok_r(p, " ,\t\n", &p))) {
				try {
					tmpVal = stof(tok); 
				}
				catch(...) {
					tmpVal = NAcode;
				}
				tmpInd.PGS.push_back(tmpVal);
			}
			if (tmpInd.PGS.size() == nTrait) {
				tmpInd.StandRes = tmpStandRes;
				SampleList[tmpIndIdx] = tmpInd;
				i++;
			}
			else {
				printf("Individual %s does not have %d PGS input, discarded.\n", tmpIndIdx.c_str(), nTrait);
			}
           tmpInd = {};
		}
	}
	printf("Read genetic risk file, done. %ld individuals has genetic risk input.\n", i);
	fclose(GenoFile);
}


void ReadPhenoFile(const char *FileName) {
	char *tok; char *p;
	long int k;
	string tmpIndIdx;
	double tmpVal;
	vector<int> PhenoIdx(nTrait, -1);
	vector<double> tmpPhenoVal(nTrait, NAcode);

	FILE *PhenoFile;
	PhenoFile = fopen(FileName, "r");

	if (PhenoFile == NULL) {
		printf("Cannot open phenotype file %s.\n", FileName);
		exit(0);
	}
	else {
		fgets(buffer, sizeof(buffer), PhenoFile); // header line
		p = buffer;
		tok = strtok_r(p, " ,\t\n", &p);
		k = 0;
		while ((tok = strtok_r(p, " ,\t\n", &p))) {
			if (TraitList.find(tok) != TraitList.end()) {
				PhenoIdx[k] = TraitList[tok].Idx;
				TraitList[tok].HasPheno = 1;
				k++;
			}
		}
		printf("Considering %d traits having both genetic risk and phenotype.\n", k);
		for (auto iter = TraitList.begin(); iter != TraitList.end();) {
			if (iter->second.HasPheno != 1)
				iter = TraitList.erase(iter);
			else
				iter++;		
		}

		while (fgets(buffer, sizeof(buffer), PhenoFile) != NULL) {
			p = buffer;
			tok = strtok_r(p, " ,\t\n", &p); //ID
			tmpIndIdx = tok;
			if (SampleList.find(tmpIndIdx) != SampleList.end()) {
				k = 0;
				while ((tok = strtok_r(p, " ,\t\n", &p))) {
					if (PhenoIdx[k] != -1) {
						try {
							tmpPhenoVal[PhenoIdx[k]] = stof(tok); 
						}
						catch(...) {
							//
						}
					}
					k++;
				}
				SampleList[tmpIndIdx].Pheno = tmpPhenoVal;
				std::fill(tmpPhenoVal.begin(), tmpPhenoVal.end(), NAcode);
			}
		}
	}
	for (auto iter = SampleList.begin(); iter != SampleList.end();) {
		if (SampleList[iter->first].Pheno.empty()) {
			iter = SampleList.erase(iter);
		}
		else {
			iter++;
		}
	}
	nSample = SampleList.size();
	printf("Read phenotype file, done. %d individuals has both genetic risk and phenotype input.\n", nSample);
	fclose(PhenoFile);
}


void ReadCovarFile(const char *FileName) {
	char *tok; char *p;
	int i = 0;
	long int n;
	string tmpIndIdx;
	vector<double> tmpCovarVal;

	FILE *CovarFile;
	CovarFile = fopen(FileName, "r");

	if (CovarFile == NULL) {
		printf("Cannot open genetic risk file %s.\n", FileName);
		exit(0);
	}
	else {
		fgets(buffer, sizeof(buffer), CovarFile); // header line
		p = buffer;
		tok = strtok_r(p, " ,\t\n", &p); // ID
		i = 0;
		while ((tok = strtok_r(p, " ,\t\n", &p))) {
			i++;
		}
		nCovar = i;

		while (fgets(buffer, sizeof(buffer), CovarFile) != NULL) {
			p = buffer;
			tok = strtok_r(p, " ,\t\n", &p); //ID
			tmpIndIdx = tok;
			if (SampleList.find(tmpIndIdx) != SampleList.end()) {
				i = 0;
				while ((tok = strtok_r(p, " ,\t\n", &p))) {
					try {
						tmpCovarVal.push_back(stof(tok)); 
					}
					catch(...) {
						break;
					}
					i++;
				}
				if (i == nCovar) {
					SampleList[tmpIndIdx].Covar = tmpCovarVal;
					tmpCovarVal.clear();
				}
				else {
					SampleList.erase(tmpIndIdx);
					tmpCovarVal.clear();
					printf("Individual %s does not have %d covariates, discarded.\n", tmpIndIdx.c_str(), nCovar);
				}
			}
		}
	}
	nSample = SampleList.size();
	printf("Read covariate file, done. Leaving %ld valid individuals each with %d covariates.\n", nSample, nCovar);
	fclose(CovarFile);
}


void MakeMat() {
	long int n;
	int i;
	GeneMat.set_size(nSample, nTrait);
	PhenoMat.set_size(nSample, nTrait);
	CovarMat.set_size(nSample, nCovar + 1);
	
	n = 0;
	for (auto iter = SampleList.begin(); iter != SampleList.end(); iter++) {
		for (i = 0; i < nTrait; i++) {
			GeneMat(n, i) = iter->second.PGS[i];
			PhenoMat(n, i) = iter->second.Pheno[i];
		}
		for (i = 0; i < nCovar; i++) {
			CovarMat(n, i) = iter->second.Covar[i];
		}
		CovarMat(n, nCovar) = 1;
		n++;
	}
	printf("Making input matrices, done.\n");
}


void GetStandRes(string Trait) {
	long int n;
	int iTrait = TraitList[Trait].Idx;
	double meanX, meanY, sdX, sdY, meanRes, sdRes;
	arma::mat X;
	arma::colvec Y;

	X = join_rows(GeneMat.col(iTrait), CovarMat);
	Y = PhenoMat.col(iTrait);

	arma::uvec AllIdx = arma::regspace<arma::uvec>(0, nSample - 1);
	arma::uvec idx_Gene = find(GeneMat.col(iTrait) == NAcode);
	arma::uvec idx_Pheno = find(PhenoMat.col(iTrait) == NAcode);
	arma::uvec idx_rm = unique(join_cols(idx_Gene,idx_Pheno));
	std::vector<int> KeepIdx;
	std::set_difference(AllIdx.begin(), AllIdx.end(), idx_rm.begin(), idx_rm.end(), std::inserter(KeepIdx, KeepIdx.begin()));

	X = X.rows(arma::conv_to<arma::uvec>::from(KeepIdx));
	Y = Y.rows(arma::conv_to<arma::uvec>::from(KeepIdx));

	long int N = X.n_rows, k = X.n_cols;
	arma::colvec coef(k);

	sdX = arma::stddev(X.col(0));
	sdY = arma::stddev(Y);

   	if (sdX != 0 && sdY != 0) {
   		// OutStatFile = fopen(OutStat,"a");
		meanX = arma::mean(X.col(0));
		meanY = arma::mean(Y);
		TraitList[Trait].Corr = arma::norm_dot( (X.col(0)-meanX)/sdX, (Y-meanY)/sdY);

		if (arma::solve(coef, X, Y, arma::solve_opts::allow_ugly + arma::solve_opts::fast)) {
			arma::colvec resid = Y - X*coef;

			double TSS = pow(arma::norm(Y, 2),2);
			double ESS = pow(arma::norm(resid, 2),2);
			TraitList[Trait].R2 = 1 - ESS/TSS;
			TraitList[Trait].adjR2 = 1 - (1 - TraitList[Trait].R2)*(N-1)/(N-k);

			meanRes = arma::mean(resid);
			sdRes = arma::stddev(resid);
			resid = (resid-meanRes)/sdRes;

			n = 0;
			for (auto iter = SampleList.begin(); iter != SampleList.end(); iter++) {
				if (SampleList[iter->first].PGS[iTrait] != NAcode && SampleList[iter->first].Pheno[iTrait] != NAcode) {
					SampleList[iter->first].StandRes[iTrait] = resid(n);
					n++;
				}
			}

			printf("Trait: %s, done.\n", Trait.c_str());
			printf("Regression Coeffs:");
			int i;
			for (i = 0; i < k; i++) {
				printf(" %lf", coef(i));
			}
			printf("\n");
			printf("#Regression Ind = %ld, geno~Pheno Corr = %lf, R2 = %lf, adjR2 = %lf\n",  
				N, TraitList[Trait].Corr, TraitList[Trait].R2, TraitList[Trait].adjR2);
		}
		else
			printf("Trait %s linear regression failed.\n", Trait.c_str());
	}
	else {
		printf("Trait %s does not have non-zero variance in both genetic risk and phenotype.\n", Trait.c_str());
	}
}


void GetAggRes() {
	double TotW;
	double wSumRes;

	for (auto iter = TraitList.begin(); iter != TraitList.end();) {
		if (iter->second.adjR2 < R2thres)
			iter = TraitList.erase(iter);
		else
			iter++;
	}

	for (auto iter1 = SampleList.begin(); iter1 != SampleList.end(); iter1++) {
		TotW = 0.0;
		wSumRes = 0.0;
		for (auto iter2 = TraitList.begin(); iter2 != TraitList.end(); iter2++) {
			if (iter1->second.StandRes[iter2->second.Idx] != NAcode) {
				wSumRes += abs(iter1->second.StandRes[iter2->second.Idx]) * iter2->second.adjR2;
				TotW += iter2->second.adjR2;
			}
		}
		SampleList[iter1->first].AggStandRes = wSumRes/TotW;
	}
}


int main(int argc, char const *argv[]) {
	nTrait = 0;
	nCovar = 0;
	nSample = 0;

	ReadParam(argv[1]);
	ReadGenoFile(InGenoFile.c_str());
	ReadPhenoFile(InPhenoFile.c_str());
	if (UseCovar)
		ReadCovarFile(InCovarFile.c_str());
	MakeMat();

	for (auto iter = TraitList.begin(); iter != TraitList.end();iter++) {
		GetStandRes(iter->first);
	}
	GetAggRes();

	FILE *Out;
	Out = fopen(OutFile.c_str(),"w");
	// Set default values
	if (Out == NULL) {
	    printf("Cannot open output file.\n");
	    exit(0);
	}
	fprintf(Out, "IndID\tAggAbsResid");
	for (auto iter = TraitList.begin(); iter != TraitList.end();iter++) {
		fprintf(Out, "\tStandResid_%s", iter->first.c_str());
	}
	fprintf(Out, "\n");

	for (auto iter1 = SampleList.begin(); iter1 != SampleList.end(); iter1++) {
		if (SampleList[iter1->first].AggStandRes != NAcode)
			fprintf(Out, "%s\t%lf", iter1->first.c_str(), SampleList[iter1->first].AggStandRes);
		else
			fprintf(Out, "%s\tNA", iter1->first.c_str());

		for (auto iter2 = TraitList.begin(); iter2 != TraitList.end();iter2++) {
			if (SampleList[iter1->first].StandRes[iter2->second.Idx] != NAcode)
				fprintf(Out, "\t%lf", SampleList[iter1->first].StandRes[iter2->second.Idx]);
			else
				fprintf(Out, "\tNA");
		}
		fprintf(Out, "\n");
	}
	fclose(Out);
}

