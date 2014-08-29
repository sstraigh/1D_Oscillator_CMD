#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#define COL_NUMBER 11

int main(){

    std::vector<double> locations;

    while (std::cin){

	for (int i=0; i<COL_NUMBER; ++i){

	    std::string holder;
	    std::cin>>holder;
	    if (i==5){
		std::istringstream iss(holder);
		double temp=0.0;
		iss>>temp;
		locations.push_back(temp);

	    }
	}

    }

    double min=0.0;
    double max=0.0;

    for (std::vector<double>::iterator it=locations.begin(); it!=locations.end(); ++it){

	if ((*it)<min) min=(*it);
	if ((*it)>max) max=(*it);

    }

    double range=max-min;
    int num_entries=(range*10);
    double bin_size=range/((double)num_entries);

    int * frequency_of_locations=new int[num_entries];
    std::fill(frequency_of_locations, frequency_of_locations+num_entries, 0);


    for (std::vector<double>::iterator it=locations.begin(); it!=locations.end(); ++it){

	for (int i=1; i<=num_entries; ++i){

	    double lower_bound=(min + ((i-1) * bin_size));
	    double upper_bound=(min + ( i * bin_size));

	    if ((*it)>lower_bound && (*it)<upper_bound) ++frequency_of_locations[i-1];

	}
    }

    for (int i=0; i<num_entries; ++i){

	std::cout<<std::setw(18)<<(min+(i*bin_size))<<" "
	    <<std::setw(18)<<frequency_of_locations[i]
	    <<'\n';
    }

    return 0;
}
