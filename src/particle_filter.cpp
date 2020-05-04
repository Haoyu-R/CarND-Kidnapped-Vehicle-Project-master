///**
// * particle_filter.cpp
// *
// * Created on: Dec 12, 2016
// * Author: Tiffany Huang
// */
//
#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

static std::default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	/**
	 * TODO: Set the number of particles. Initialize all particles to
	 *   first position (based on estimates of x, y, theta and their uncertainties
	 *   from GPS) and all weights to 1.
	 * TODO: Add random Gaussian noise to each particle.
	 * NOTE: Consult particle_filter.h for more information about this method
	 *   (and others in this file).
	 */

	std::normal_distribution<double> gps_x(x, std[0]);
	std::normal_distribution<double> gps_y(y, std[1]);
	std::normal_distribution<double> theta_(theta, std[2]);

	num_particles = 100;  // TODO: Set the number of particles

	for (int i = 0; i < num_particles;i++) {
		Particle particle;
		particle.id = i;
		particle.x = gps_x(gen);
		particle.y = gps_y(gen);
		particle.theta = theta_(gen);
		particle.weight = 1;
		particles.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
	double velocity, double yaw_rate) {
	/**
	 * TODO: Add measurements to each particle and add random Gaussian noise.
	 * NOTE: When adding noise you may find std::normal_distribution
	 *   and std::default_random_engine useful.
	 *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	 *  http://www.cplusplus.com/reference/random/default_random_engine/
	 */
	std::normal_distribution<double> x(0, std_pos[0]);
	std::normal_distribution<double> y(0, std_pos[1]);
	std::normal_distribution<double> theta(0, std_pos[2]);

	for (auto i = particles.begin(); i != particles.end();i++)
	{
		if (yaw_rate < 0.0001) {
			i->x = i->x + velocity * cos(i->theta) * delta_t + x(gen);
			i->y = i->y + velocity * sin(i->theta) * delta_t + y(gen);
		}
		else {
			i->x = i->x + velocity / yaw_rate * (sin(i->theta + yaw_rate * delta_t) - sin(i->theta)) + x(gen);
			i->y = i->y + velocity / yaw_rate * (-cos(i->theta + yaw_rate * delta_t) + cos(i->theta)) + y(gen);
		}
		i->theta = i->theta + delta_t*yaw_rate + theta(gen);
	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
	vector<LandmarkObs>& observations) {
	/**
	 * TODO: Find the predicted measurement that is closest to each
	 *   observed measurement and assign the observed measurement to this
	 *   particular landmark.
	 * NOTE: this method will NOT be called by the grading code. But you will
	 *   probably find it useful to implement this method and use it as a helper
	 *   during the updateWeights phase.
	 */
	for (auto i = observations.begin();i != observations.end();i++) {
		int id = -1;
		double temp_dist = 50;
		for (auto j = predicted.begin();j != predicted.end();j++) {
			double dist_ = dist(i->x, i->y, j->x, j->y);
			if (dist_ < temp_dist)
			{
				temp_dist = dist_;
				id = j->id;
			}
		}
		i->id = id;


	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const vector<LandmarkObs>& observations,
	const Map& map_landmarks) {

	

	double weight_sum = 0;
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	for (auto i = particles.begin();i != particles.end();i++)
	{

		double p_x = i->x;
		double p_y = i->y;
		double p_theta = i->theta;


		vector<LandmarkObs> predictions;
		for (auto j = map_landmarks.landmark_list.begin();j != map_landmarks.landmark_list.end();j++)
		{
			double d = dist(j->x_f, j->y_f, i->x, i->y);
			if (d > sensor_range)
			{
				continue;
			}

			LandmarkObs lmo;
			lmo.id = j->id_i;
			lmo.x = j->x_f;
			lmo.y = j->y_f;

			predictions.push_back(lmo);
		}

		//vector<LandmarkObs> transformed_observation;
		//for (size_t j = 0; j < observations.size(); ++j) {

		//	
		//	// Convert observation from particle(vehicle) to map coordinate system
		//	LandmarkObs rototranslated_obs;
		//	rototranslated_obs.x = cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y + p_x;
		//	rototranslated_obs.y = sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y + p_y;

		//	transformed_observation.push_back(rototranslated_obs);
		//}

		vector<LandmarkObs> transformed_observation;
		for (auto z = observations.begin();z != observations.end();z++)
		{

			double old_x = z->x;
			double old_y = z->y;

			double new_x = p_x + cos(p_theta) * old_x - sin(p_theta) * old_y;
			double new_y = p_y + sin(p_theta) * old_x + cos(p_theta) * old_y;

			LandmarkObs lmo;
			lmo.x = new_x;
			lmo.y = new_y;
			transformed_observation.push_back(lmo);
		}

		dataAssociation(predictions, transformed_observation);

	//	double particle_likelihood = 1.0;

	//	double mu_x, mu_y;
	//	for (const auto& obs : transformed_observation) {

	//		// Find corresponding landmark on map for centering gaussian distribution
	//		for (const auto& land : predictions)
	//			if (obs.id == land.id) {
	//				mu_x = land.x;
	//				mu_y = land.y;
	//				break;
	//			}

	//		double norm_factor = 2 * M_PI * std_x * std_y;
	//		double prob = exp(-(pow(obs.x - mu_x, 2) / (2 * std_x * std_x) + pow(obs.y - mu_y, 2) / (2 * std_y * std_y)));

	//		particle_likelihood *= prob / norm_factor;
	//	}

	//	i->weight = particle_likelihood;

	//} // end loop for each particle

	//// Compute weight normalization factor
	//double norm_factor = 0.0;
	//for (const auto& particle : particles)
	//	norm_factor += particle.weight;

	//// Normalize weights s.t. they sum to one
	//int iter = 0;
	//	for (auto& particle : particles)
	//{
	//	particle.weight /= (norm_factor);
	//	weights[iter] = particle.weight;
	//	iter++;
	//}

		for (auto j = transformed_observation.begin(); j != transformed_observation.end();j++)
		{
			double mu_x;
			double mu_y;


			i->weight *= Gaussian2D(j->x, j->y, (map_landmarks.landmark_list[j->id-1]).x_f, (map_landmarks.landmark_list[j->id-1]).y_f, std_landmark[0], std_landmark[1]);
			//i->weight *= Gaussian2D(j->x, j->y, mu_x, mu_y, std_landmark[0], std_landmark[1]);
		
		}

		weight_sum += i->weight;
	}

	int iter = 0;
	for (auto j = particles.begin();j != particles.end();j++)
	{
		j->weight /= weight_sum;
		weights[iter] = j->weight;
		iter++;
	}
}

//void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
//	const vector<LandmarkObs>& observations,
//	const Map& map_landmarks) {
//
//	// Gather std values for readability
//	double std_x = std_landmark[0];
//	double std_y = std_landmark[1];
//
//	// Iterate over all particles
//	for (size_t i = 0; i < num_particles; ++i) {
//
//		// Gather current particle values
//		double p_x = particles[i].x;
//		double p_y = particles[i].y;
//		double p_theta = particles[i].theta;
//
//		// List all landmarks within sensor range
//		vector<LandmarkObs> predicted_landmarks;
//
//		for (const auto& map_landmark : map_landmarks.landmark_list) {
//			int l_id = map_landmark.id_i;
//			double l_x = (double)map_landmark.x_f;
//			double l_y = (double)map_landmark.y_f;
//
//			double d = dist(p_x, p_y, l_x, l_y);
//			if (d < sensor_range) {
//				LandmarkObs l_pred;
//				l_pred.id = l_id;
//				l_pred.x = l_x;
//				l_pred.y = l_y;
//				predicted_landmarks.push_back(l_pred);
//			}
//		}
//
//		// List all observations in map coordinates
//		vector<LandmarkObs> observed_landmarks_map_ref;
//		for (size_t j = 0; j < observations.size(); ++j) {
//
//			// Convert observation from particle(vehicle) to map coordinate system
//			LandmarkObs rototranslated_obs;
//			rototranslated_obs.x = cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y + p_x;
//			rototranslated_obs.y = sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y + p_y;
//
//			observed_landmarks_map_ref.push_back(rototranslated_obs);
//		}
//
//		// Find which observations correspond to which landmarks (associate ids)
//		dataAssociation(predicted_landmarks, observed_landmarks_map_ref);
//
//		// Compute the likelihood for each particle, that is the probablity of obtaining
//		// current observations being in state (particle_x, particle_y, particle_theta)
//		double particle_likelihood = 1.0;
//
//		double mu_x, mu_y;
//		for (const auto& obs : observed_landmarks_map_ref) {
//
//			// Find corresponding landmark on map for centering gaussian distribution
//			for (const auto& land : predicted_landmarks)
//				if (obs.id == land.id) {
//					mu_x = land.x;
//					mu_y = land.y;
//					break;
//				}
//
//			double norm_factor = 2 * M_PI * std_x * std_y;
//			double prob = exp(-(pow(obs.x - mu_x, 2) / (2 * std_x * std_x) + pow(obs.y - mu_y, 2) / (2 * std_y * std_y)));
//
//			particle_likelihood *= prob / norm_factor;
//		}
//
//		particles[i].weight = particle_likelihood;
//
//	} // end loop for each particle
//
//	// Compute weight normalization factor
//	double norm_factor = 0.0;
//	for (const auto& particle : particles)
//		norm_factor += particle.weight;
//
//	// Normalize weights s.t. they sum to one
//	int iter = 0;
//	for (auto& particle : particles)
//	{
//		particle.weight /= (norm_factor);
//		weights[iter] = particle.weight;
//		iter++;
//	}
//
//	
//}


void ParticleFilter::resample() {

	std::discrete_distribution<int> d(weights.begin(), weights.end());

	vector<Particle> new_particles;
	for (int i = 0;i < weights.size();i++)
	{
		int number = d(gen);
		particles[i].weight = 1;
		new_particles.push_back(particles[number]);
	}

	particles = new_particles;
	
}


void ParticleFilter::SetAssociations(Particle& particle,
	const vector<int>& associations,
	const vector<double>& sense_x,
	const vector<double>& sense_y) {
	// particle: the particle to which assign each listed association, 
	//   and association's (x,y) world coordinates mapping
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates
	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
	vector<double> v;

	if (coord == "X") {
		v = best.sense_x;
	}
	else {
		v = best.sense_y;
	}

	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}


