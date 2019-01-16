/*                                                                         80->|
 * particle_filter.cpp
 *
 *  Created on: May 15, 2017
 *      Author: James William Dunn
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
using namespace std;

const int NUM_PARTICLES = 50; // min:4, max:900

// The following can be set to 1.0 or 1/n
const double INIT_WEIGHT = 1.0 / NUM_PARTICLES; 

/*
 * Approximate the distance between two 2D points within 4% error.
 * @param (x1,y1) x and y coordinates of first point
 * @param (x2,y2) x and y coordinates of second point
 * @output Approximate distance between the points
 * Performance: 2.336 secs when n=50
 */
inline double distAB(double x1, double y1, double x2, double y2) {
  double a = abs(x2-x1);
  double b = abs(y2-y1);
  // Implements alpha max beta min algorithm
  if (a>=b) return 0.960448467 * a + 0.397847893 * b;
  else return 0.960448467 * b + 0.397847893 * a;
}
/*
 * Compute the square distance between two 2D points.
 * @param (x1,y1) x and y coordinates of first point
 * @param (x2,y2) x and y coordinates of second point
 * @output The square of the distance between the points
 * Performance: 1.908 secs when n=50
 */
inline double distSqr(double x1, double y1, double x2, double y2) {
  double a = x2-x1;
  double b = y2-y1;
  return a * a + b * b;
}
/*
 * Initialize a particle filter object.
 * @param (x,y) approximate x and y coordinates of vehicle
 * @param theta Yaw of the vehicle
 * @param std Standard deviation of the measurements
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles.
  num_particles = NUM_PARTICLES;
  particles.resize(num_particles);
  
  // Create a pseudo-random number object
	default_random_engine rnd;
  /*cout << "x:" << x << ", y:" << y << ", theta:" << theta 
       << ", std_x:" << std[0] << ", std_y:" << std[1] 
       << ", std_theta:" << std[2] << ", rnd():" << rnd() << endl;
       */
  
  // Prepare to add random Gaussian noise to each particle.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_t(theta, std[2]);
  
	// Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) 
  // and all weights to 1.0/num_particles.
  for (int i = 0; i < num_particles; i++) {
		particles[i].id = i;
		particles[i].x = dist_x(rnd);
		particles[i].y = dist_y(rnd);
		particles[i].theta = dist_t(rnd);
		particles[i].weight = INIT_WEIGHT;
    // cout << "  particle is here:" << particles[i].x << "\t" 
    // << particles[i].y << "\t" << particles[i].theta << endl;
	}
  is_initialized = true;
}

/*
 * Predict the vehicle's next state.
 * @param delta_t in seconds (e.g. 0.1)
 * @param std_pos Standard deviation of the measurements
 * @param velocity in meters/sec
 * @param yaw_rate in radians/sec
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], 
       double velocity, double yaw_rate) {
	// Add measurements to each particle
	// en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// cplusplus.com/reference/random/default_random_engine
  // cout << "prediction:" << delta_t << ":" << std_pos[0] << ":" 
  // << velocity << ":" << yaw_rate << endl;
  default_random_engine rnd;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_t(0, std_pos[2]);

  for (int i = 0; i < num_particles; i++) {
    // Avoid division by zero
    if (fabs(yaw_rate) < 1e-5) {
      particles[i].x += delta_t * velocity * cos(particles[i].theta);
      particles[i].y += delta_t * velocity * sin(particles[i].theta);
    } else {
      particles[i].x += velocity / yaw_rate * 
        (sin(particles[i].theta + yaw_rate * delta_t) -
        sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * 
        (-cos(particles[i].theta + yaw_rate * delta_t) + 
        cos(particles[i].theta));
    }
    // Add random Gaussian noise.
    particles[i].x += dist_x(rnd);
    particles[i].y += dist_y(rnd);
    particles[i].theta += delta_t * yaw_rate + dist_t(rnd);
    // cout << "  particle is here:" << particles[i].x << "\t" 
    // << particles[i].y << endl;
  }
}

// Unused function: method is implemented inline below
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
                std::vector<LandmarkObs>& observations) {
	// Find the landmark that is closest to each observed measurement
}

/*
 * Update the particle weights
 * @param sensor_range in meters (e.g. 50)
 * @param std_landmark Standard deviation of the landmark
 * @param observations in the VEHICLE'S coordinate system
 * @param map_landmarks in the MAP'S coordinate system
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian 
  // distribution. 
	//   en.wikipedia.org/wiki/Multivariate_normal_distribution
	//   willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   planning.cs.uiuc.edu/node99.html
  
  //cout << "updateWeights..." << endl;
  if (observations.size() == 0) return; // no observations, bail out
  
  double sig_denom = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
  double sig_x2 = std_landmark[0] * std_landmark[0];
  double sig_y2 = std_landmark[1] * std_landmark[1];
  
  // Update all particles
  for (int i = 0; i < particles.size(); i++) {
  
    auto transformed = observations;
    // Walk through the observations list
    for (int j = 0; j < observations.size(); j++) {
      // Convert vehicle coords to world coords: rotation then translation
      // cout << "Observed:" << observations[j].x << "\t" << observations[j].y
      // << endl;
      double cos_theta = cos(particles[i].theta);
      double sin_theta = sin(particles[i].theta);
      transformed[j].x = observations[j].x * cos_theta - 
        observations[j].y * sin_theta + particles[i].x;
      transformed[j].y = observations[j].x * sin_theta + 
        observations[j].y * cos_theta + particles[i].y;
      
      double lowest = sensor_range;
      // Compare observations to landmarks within sensor range
      for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
        
        double d_AB = distSqr(transformed[j].x, transformed[j].y, 
          map_landmarks.landmark_list[k].x_f, 
          map_landmarks.landmark_list[k].y_f);
        //if (d_AB < sensor_range)
        //  cout << "  landmark:" << k+1 << ", distance:"<< d_AB << endl;
        if (d_AB < lowest) { // find closest
          lowest = d_AB;
          // landmark ids range from 1 to 42
          transformed[j].id = map_landmarks.landmark_list[k].id_i-1; 
        }
      }
      //cout << "  landmark:" << transformed[j].id+1 << ", closest:" 
      //<< lowest << endl;
      
      // Compute a new weight based on distance to observed landmarks
      double x = pow(transformed[j].x - 
        map_landmarks.landmark_list[transformed[j].id].x_f, 2) / (2 * sig_x2);
      double y = pow(transformed[j].y - 
        map_landmarks.landmark_list[transformed[j].id].y_f, 2) / (2 * sig_y2);
      double weight = exp(-(x + y)) * sig_denom;
      
      // Multiply all the weights to compute total probability of particle
      if (weight > 0.0) particles[i].weight *= weight;
      //cout << "  id:" << particles[i].id << ", w:" << particles[i].weight 
      //<< endl;
    }
  }
  //cout << "  ...updateWeights" << endl;
}

/*
 * Resample the particles
 */
void ParticleFilter::resample() {
	// Replace with probability proportional to particle weight. 
	// en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  // cout << "resample..." << endl;
  
  // output prior to resampling for alternate visualization
  //write("data/particles.txt");
  vector<Particle> new_p;
  default_random_engine rnd;
	// Put all weights into a vector
	for (int i = 0; i < num_particles; i++)
		weights.push_back(particles[i].weight);
  
  discrete_distribution<int> d_d(weights.begin(), weights.end());
  // Resample particles with probability proportional to weight
  for (int i = 0; i < num_particles; i++) {
    int index = d_d(rnd);
    new_p.push_back(particles[index]);
  }
  
  // Replace particles and reset weights
  // could do this: particles = new_p;
  // however, also need to reset weight
  for (int i = 0; i < num_particles; i++) {
    particles[i].x = new_p[i].x;
    particles[i].y = new_p[i].y;
    particles[i].theta = new_p[i].theta;
    particles[i].weight = INIT_WEIGHT;  // set up for the next cycle
    // cout << "  particle is here:" << particles[i].x << "\t" 
    // << particles[i].y << "\t" << particles[i].theta << endl;
    }
  weights.resize(0);
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " 
      << particles[i].theta << "\n";
	}
	dataFile.close();
}

/* The following is for the Udacity Particle Filter Project Visualizer

See the video for output from the test run: https://vimeo.com/217426713

#include <sstream>
#include <iterator>
Particle ParticleFilter::SetAssociations(Particle particle, 
  std::vector<int> associations, 
  std::vector<double> sense_x, std::vector<double> sense_y) {
	//particle: the particle to assign each listed association, and association's
  // (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best) {
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best) {
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


*/