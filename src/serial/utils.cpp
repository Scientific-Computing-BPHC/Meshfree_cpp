#include "utils.hpp"

std::string readFile(char* file_name)
{
	std::string gridfile, tp;

	std::fstream datafile(file_name, ios::in);
	if(datafile.is_open())
	{
		cout<<"File opened"<<endl;
		//read the file
		while(getline(datafile, tp))
		{
			//cout<<"Data File: \n"<<tp<<endl;
			gridfile.append(tp);
			gridfile.append("\n");
		}
	}
	//cout<<"Grid file: \n"<<gridfile<<endl;
	datafile.close();

	return gridfile;
}

void checkFileRead(std::vector<vec_doub> result_doub, int numPoints)
{
	std::ofstream out("checkFileRead.txt");
	for(int j=0; j<numPoints; j++)
    {
		for (int i=0; i<result_doub[j].size(); i++)
			out<<std::fixed<<std::setprecision(20)<<result_doub[j][i]<<"   ";
		out<<endl;
	}
	out.close();

}


void printPrimal(double primal[4])
{
	cout<<"Printing Primal: "<<endl;
	for(int i=0; i<4; i++)
		cout<<std::fixed<<std::setprecision(20)<<primal[i]<<" ";
	cout<<"\n";
}


