/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

using namespace std;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 500;

	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	random_device rd;
	default_random_engine gen(rd());

	// This line creates a normal (Gaussian) distribution for x, y, theta.
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; i++) {

		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);                //gt_data plus noise for each particle
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
		weights.push_back(particle.weight);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	random_device rd;
	default_random_engine gen(rd());

	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	for (int i = 0; i < num_particles; i++) {
		Particle &particle = particles[i];
		double delta_theta = yaw_rate * delta_t;

		particle.x = particle.x + (velocity / yaw_rate) * (sin(particle.theta + delta_theta) - sin(particle.theta));
		particle.y = particle.y + (velocity / yaw_rate) *(cos(particle.theta) - cos(particle.theta + delta_theta));

		particle.theta = particle.theta + delta_theta;

		//adding noise to the particle predction outcome 
		normal_distribution<double> dist_x(particle.x, std_x);
		normal_distribution<double> dist_y(particle.y, std_y);
		normal_distribution<double> dist_theta(particle.theta, std_theta);

		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
	}

}

/*
 * This function is not used, instead non member function data_association_per_particle is used to perform data association
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  for (int i = 0; i<num_particles; i++)
  {
	Particle& particle = particles[i];
	vector<LandmarkObs> transformed_landmarks;     
	std::vector<LandmarkObs> predicteds_observations;

    double px = particle.x;
    double py = particle.y;
    double ptheta = particle.theta;

	transform_landmarks_per_particle(particle, map_landmarks, transformed_landmarks);

    //std::sort(transformed_landmarks.begin(), transformed_landmarks.end());

	data_association_per_particle(transformed_landmarks, observations, particle, predicteds_observations);

	particle.weight = measurement_prob(predicteds_observations, std_landmark, observations);
	weights[i] = particle.weight;
  } // End particle loop
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::vector<Particle> new_particles;

  double max_weight = *max_element(begin(weights), end(weights));

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> ddouble(0.0, 2.0*max_weight);
  std::uniform_int_distribution<int> dint(0, num_particles);

  // draw with some noise
  std::normal_distribution<double> dx(0,0.003);
  std::normal_distribution<double> dy(0,0.003);
  std::normal_distribution<double> dtheta(0,0.0001);

  int idx = dint(gen);
  double beta = 0.0;

  for (int i = 1; i<num_particles; i++)
  {
    beta = beta + ddouble(gen);
    while (beta > particles[idx].weight)
    {
      beta = beta - particles[idx].weight;
      idx = (idx+1)%num_particles;
    }
    Particle particle;
    particle.x = particles[idx].x+dx(gen);
    particle.y = particles[idx].y+dy(gen);
    particle.theta = particles[idx].theta+dtheta(gen);
    //particle.x = particles[idx].x;
    //particle.y = particles[idx].y;
    //particle.theta = particles[idx].theta;

    particle.weight = particles[idx].weight;
    new_particles.push_back(particle);
  }
  particles = new_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}


/**
* Transform the map landmarks to the particle coordinate system.
* @param particle The particle whose coordinate system defines the transformation
* @param map_landmarks Map class containing map landmarks
* @output transformed_landmarks， the vector of landmarks transformed to the particle coordinate system
*/
void ParticleFilter::transform_landmarks_per_particle(const Particle particle, const Map& map_landmarks, vector<LandmarkObs>& transformed_landmarks) {

	for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
		const Map::single_landmark_s& single_landmark = map_landmarks.landmark_list[i];
		LandmarkObs transformed_landmark;
		transformed_landmark.id = single_landmark.id_i;
		double cos_theta = cos(particle.theta - M_PI / 2);
		double sin_theta = sin(particle.theta - M_PI / 2);
		transformed_landmark.x = -(single_landmark.x_f - particle.x) * sin_theta + (single_landmark.y_f - particle.y) * cos_theta;
		transformed_landmark.y = -(single_landmark.x_f - particle.x) * cos_theta - (single_landmark.y_f - particle.y) * sin_theta;
		transformed_landmarks.push_back(transformed_landmark);

	}
	return;
}

/**
 * Associate each observation to its mostly likely predicted landmark measurements for each particle
 * @param transformed_landmarks, the vector of landmarks transformed to the particle coordinate system
 * @param observations, the list of actual measurements
 * @param particle, the particle being processed
 * @output predicteds_observations, the vector of predicted measurements
 */

void ParticleFilter::data_association_per_particle(const std::vector<LandmarkObs> transformed_landmarks,
	const std::vector<LandmarkObs> observations,
	const Particle particle,
	std::vector<LandmarkObs>& predicteds_observations) {


	for (int i = 0; i < observations.size(); i++) {
		const LandmarkObs& obs = observations[i];

		//Finding the closest landmark for the transformed observation 
		//and add it as the associated landmark.
		double clostest_landmark_dist = -1;
		int clostest_landmark_ind = -1;

		for (int j = 0; j < transformed_landmarks.size(); j++) {
			const LandmarkObs& single_landmark = transformed_landmarks[j];
			double x_dist = single_landmark.x - obs.x;
			double y_dist = single_landmark.y - obs.y;
			double dist = sqrt(x_dist*x_dist + y_dist*y_dist);

			if (clostest_landmark_dist == -1 || dist < clostest_landmark_dist) {
				clostest_landmark_dist = dist;
				clostest_landmark_ind = j;
			}
		}
		const LandmarkObs& predicteds_obs = transformed_landmarks[clostest_landmark_ind];
		predicteds_observations.push_back(predicteds_obs);
	}

}

//calculates how likely a measurement should be
double ParticleFilter::measurement_prob(const std::vector<LandmarkObs>& predicteds_observations,
	double std_landmark[], const std::vector<LandmarkObs>& observations) {
	// Reference: Udacity SDC, Term 2 Lesson 13- Particle Filters. Section 14 - Exercise Importance Weight

	long double prob = 1.0;
	for (int i = 0; i < observations.size(); i++) {

		const LandmarkObs& measurement_obs = observations[i];
		const LandmarkObs& predicted_obs = predicteds_observations[i];

		//multiply density for all predicted measurements
		prob *= Gaussian(predicted_obs, std_landmark, measurement_obs);

	}
	return prob;
}

//# calculates the probability of x,y for 2-dim Gaussian with mean mu and var. sigma
double ParticleFilter::Gaussian(const LandmarkObs predicted_obs,
	double std_landmark[],
	const LandmarkObs measurement_obs){
	// Reference: Udacity SDC, Term 2 Lesson 13- Particle Filters. Section 14 - Exercise Importance Weight	

	double x = measurement_obs.x;
	double y = measurement_obs.y;

	double mu_x = predicted_obs.x;
	double mu_y = predicted_obs.y;

	double gap_x = x - mu_x;
	double gap_y = y - mu_y;

	long double denomin = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
	long double x_nom = ((x - mu_x) * (x - mu_x)) / (std_landmark[0] * std_landmark[0]);
	long double y_nom = ((y - mu_y) * (y - mu_y)) / (std_landmark[1] * std_landmark[1]);
	long double expon = exp(-0.5*(x_nom + y_nom));
	long double prob = denomin * expon;

	return prob;
}
