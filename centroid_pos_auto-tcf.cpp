#include <iomanip>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#define SQUARE yes
//#define CUBE yes


//usage: files are read in from cin at the command line using '<' operator,
//       output files are labeled similarly using '>'

int main(int argc, char * argv[]){

     if (argc<4){
       	 std::cerr<<"Usage: Create-Position-Autocorrelation column_containing_positions total_number_of_columns number_of_data_points < input.dat > output.dat"<<std::endl;
	 return 1;
     }

     std::istringstream ss(argv[1]);
     int column_containing_positions=0;
     ss>>column_containing_positions;

     std::istringstream sss(argv[2]);
     int total_column_number=0;
     sss>>total_column_number;

     std::istringstream ssss(argv[3]);
     int num_data=0;
     ssss>>num_data;

     double * positions=new double[num_data];
     double * time = new double[num_data];
 
     std::fill(positions, positions+num_data, 0.0);
     std::fill(time, time+num_data, 0.0);

     int count=0;

     while (std::cin){

	 for (int i=1; i<=total_column_number; ++i){

	     std::string holder;
	     std::cin>>holder;

	     if (i==1){
		 double time_temp=0.0;
		 std::istringstream iss(holder);
		 iss>>time_temp;
		 time[count]=time_temp;
	     }

	     else if (i==column_containing_positions){
		 double temp=0.0;
		 std::istringstream iss(holder);
		 iss>>temp;
		 positions[count]=temp;
	     }
	 }

	 ++count;

     }

     double * corr_func=new double[num_data];
     std::fill(corr_func, corr_func+num_data, 0.0);

     int * num_entries=new int[num_data];
     std::fill(num_entries, num_entries+num_data, 0);

     double * corr_norm=new double[num_data];
     std::fill(corr_norm, corr_norm+num_data, 0.0);

     for (int i=0; i<num_data; ++i){

	 for (int k=i; k<num_data; ++k){

	     double pos_time_0=positions[i];
	     double pos_time_t=positions[k];
	     int index=(k-i);

	    //	     if (index<0) index*=-1;

	     corr_func[index]+=(pos_time_0*pos_time_0*pos_time_t*pos_time_t);
	     ++num_entries[index];
	     corr_norm[index]+=(pos_time_0*pos_time_0* pos_time_0*pos_time_0);
	 }
     }

     double max_val=0.0;


     //enough data points have been sampled that the last of the autocorrelation
     //function can be truncated, as there are not enough averages there to let 
     //the qtcf converge

     for (int i=0; i<(num_data/5); ++i){
	 corr_func[i]/=num_entries[i];
	 if (corr_func[i]>max_val) max_val=corr_func[i];
     }
     

     //corr_func is printed out in its normalized and averaged form
     for (int i=0; i<(num_data/5); ++i){

	 std::cout<<
	     std::setw(18)<<time[i]<<
	     std::setw(18)<<(corr_func[i]/max_val)<<
	     std::setw(18)<<positions[i]<<
	     std::endl;

     }

     return 0;
}
