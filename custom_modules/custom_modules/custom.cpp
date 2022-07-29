/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

static double four_thirds_pi =  4.188790204786391;
// ellipse radius
double ellipse_radius= 1.0;
// ellipse matrix
std::vector<std::vector<double>> matrix{{1,0,0},{0,1,0},{0,0,1}};
//
std::vector<double> a_axis(3,0.0);
std::vector<double> b_axis(3,0.0);

void convert_eccentricity_to_axis(Cell* pCell)
{
	double new_volume=pCell->get_total_volume();
	double pi = 3.141592653589793238462643383279502884;
	double semimajor = parameters.doubles("major_axis_2a")/2;
	double ecc = parameters.doubles("eccentricity");
	//double vol = parameters.doubles("starting_cell_volume");
	double vol=pCell->get_total_volume();
	double b_axis_calc = semimajor*pow( (1-pow(ecc,2)), 0.5);
	//double c_axis_calc = (3*vol)/( 4*pi*pow(semimajor,2)*pow((1-ecc),0.5));
	double c_axis_calc = (3*vol)/( 4*pi*semimajor*b_axis_calc);

	/*
	std::cout << "ecc " << ecc  << std::endl;
	std::cout << "2aax " << parameters.doubles("major_axis_2a") << std::endl;
	std::cout << "aax " << semimajor  << std::endl; 
	std::cout << "bax " << b_axis_calc  << std::endl; 
	std::cout << "cax " << c_axis_calc << std::endl; 
	*/

	pCell->custom_data["axis_a"] = semimajor;
	pCell->custom_data["axis_b"] = b_axis_calc;
	pCell->custom_data["axis_c"] = c_axis_calc;
	return; 

}


void custom_update_radius( Cell* pCell, Phenotype& phenotype, double dt )
{
	double new_volume=pCell->get_total_volume();
	// 

	return; 
}

void create_cell_types( void )
{
	


	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}
double custom_volume_update(double a, double b, double c )
{
	double ellipsoid_volume= four_thirds_pi*a*b*c;
	return ellipsoid_volume;
}
void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// draw 1 cell 
	
	Cell* pC;
	
	for( int k=0; k < 1 ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < 1; n++ )
		{
      pC = create_cell( *pCD ); 
			std::vector<double> position (3,0.0);
			position[0] =pC->custom_data["cx"];
            position[1] =pC->custom_data["cy"];
            position[2] =pC->custom_data["cz"];
            //position[1] = parameters.doubles("cy"); 
			pC->assign_position( position );
			convert_eccentricity_to_axis(pC);
			pC->custom_data["rotation_about_z_axis"]=30;//in degrees
			//resize
			double new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
			//pC->set_total_volume(new_volume);
			std::cout << "vol " << new_volume<<std::endl;;
			std::cout << "aax " << pC->custom_data["axis_a"]<<std::endl;;
			std::cout << "bax " << pC->custom_data["axis_b"]<<std::endl;;
			std::cout << "cax " << pC->custom_data["axis_c"]<<std::endl;;
		}
	}
	std::cout << "test" << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }
void rotate_by_vector()
{
	return;
}
std::vector<std::vector<double>> Matrix_multiplication(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
		if(A[0].size()!=B.size())
		{
			std::cout<<" Error: No Bueno";
			return matrix;
		}

	int rows=A.size();
	int columns=B[1].size();
	std::vector<double> blank_column(columns,0.0);
	std::vector<std::vector<double>> C;
	for(int m=0; m<columns; m++)
	{
		C.push_back(blank_column);
	}
	
	/*matrix multiples A and B in the normal way, I assume you will input
	a valid matrix multiplication but ill add error checking in the future"*/
	//the usual
	//A[1][1]*B[1][1]+A[1][2]*B[2][1]+...
	//A*B
	for( int i = 0 ; i < C.size(); i++ )
	{
		for (int j = 0; j<C[1].size(); j++)
		{
			for (int k=0; k< B.size();k++)
			{
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	} 
	return C;
}


std::vector<std::vector<double>> ellipsoid_to_matrix(double axis_x,double axis_y, double axis_z)
{
	std::vector<std::vector<double>> E {{axis_x,0,0},{0,axis_y,0},{0,0,axis_z}};
	return E;
}
void transform_by_matrix( std::vector<std::vector<double>> ellipsoid,int mode=0)
{
	// mode 0, scale
	// mode 1, rotate
	// mode 2, translate
	return;
}
void neighbor_interaction(Cell *pCell_v, Cell* pCell_c)
{
	//get point(s) of intersection- if just 1 easy, if 2...
	//deepest point is furthest point into v on c that is between points of intersection
	//make vector that goes from pCell_c center through deepest point
	//calculate force and add it along that vector
	//component of force along the long axis is translation
	//component of force in other 2 cardinals do rotation
	return;

}

/*void find_longest_axis(Cell* pCell)
{

	if(pCell->custom_data["axis_a"]>= pCell->custom_data["axis_b"] &&pCell->custom_data["axis_a"] >= pCell->custom_data["axis_c"])
	{
		pCell->custom_data["longest axis"];
		return; 
	}
	if(pCell->custom_data["axis_b"]>= pCell->custom_data["axis_c"] &&pCell->custom_data["axis_b"] >= pCell->custom_data["axis_a"])
	{
		return; 
	}
	if (two-D==true)
	{
		if(pCell->custom_data["axis_c"]>= pCell->custom_data["axis_a"] &&pCell->custom_data["axis_c"] >= pCell->custom_data["axis_a"])
		{
			return; 
		}
	}
}
*/

void elongation(Cell* pCell)
{
	//find longest axis
	// longest_axis =longest_axis * scale factor

}
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	pCell->custom_data["rotation_about_z_axis"]= PhysiCell_globals.current_time*10;
	//custom_volume_update(pCell->custom_data["axis_a"], pCell->custom_data["axis_b"], pCell->custom_data["axis_c"]);

}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 




// void update_axis( Cell* pCell, Phenotype& phenotype, double dt )
// {
// 	static double four_thirds_pi =  4.188790204786391;
// 	double new_vol = phenotype.volume.total; 
// 	double scale_fac = new_vol - four_thirds_pi*pCell->custom_data["axis_a"]*pC->custom_data["axis_b"]*pC->custom_data["axis_c"]; 
// 	scale_fac = pow( scale_fac , 0.333333333333333333333333333333333333333 ); 
// 	pC->custom_data["axis_a"] *= scale_fac;
// 	pC->custom_data["axis_b"] *= scale_fac;
// 	pC->custom_data["axis_c"] *= scale_fac;
// 	return; 
// }


void custom_chemotaxis_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// bias direction is gradient for the indicated substrate 
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(phenotype.motility.chemotaxis_index);

	// change cell orientation
	custom_assign_orientation(pCell, phenotype, dt);

	// if orientation is aligned to bias
	// a.b/|a.b|
	double dot_prod = 0;
    for (int i = 0; i < 3; i++) {
   		dot_prod += pCell->state.orientation[i] * phenotype.motility.migration_bias_direction[i];
	}

	if ((dot_prod) <= 1) {
		// move up or down gradient based on this direction 
		phenotype.motility.migration_bias_direction *= phenotype.motility.chemotaxis_direction; 
		
		// normalize 
		normalize( &( phenotype.motility.migration_bias_direction ) );
	} else {
		custom_assign_orientation(pCell, phenotype, dt);
	}
	return;
}

void custom_assign_orientation(Cell* pCell, Phenotype& phenotype, double dt_)
{
	// later - set true/false for which method

	// get bias direction and turn instantaneously
	pCell->state.orientation[0] = phenotype.motility.migration_bias_direction[0];
	pCell->state.orientation[1] = phenotype.motility.migration_bias_direction[1];
	pCell->state.orientation[2] = 1;

	// get bias direction and turn incrementaly 
	pCell->state.orientation[0] += 0.1*(phenotype.motility.migration_bias_direction[0]-pCell->state.orientation[0]);
	pCell->state.orientation[1] += 0.1*(phenotype.motility.migration_bias_direction[1]-pCell->state.orientation[1]);
	pCell->state.orientation[2] = 1;

	normalize( &( pCell->state.orientation ) );

    return;
}