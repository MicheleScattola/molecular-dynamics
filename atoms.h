#ifndef ATOMS_H
#define ATOMS_H

#include <string>
#include <math.h>
#include <vector>
#include <array>
#include <stdexcept>

using std::vector;
using std::array;
using std::string;

double square(double val){ return val*val; }
double cube(double val){ return val*val*val; }

class Atom {

  public:

    // empty constructor
    Atom(){
      _mass = 1.0;
      _kSpring = 1.0;
      _cRepulsive = 0.001;
      _type = "C";
    }
    // full constructor
    Atom(string type,
      std::array<double,3> pos, std::array<double,3> vel, std::array<double,3> force,
      double m, double k, double repulsive){
        _type = type;

        // vector dimension check
        if(pos.size() !=3 || vel.size() !=3 || force.size() !=3){
          throw std::invalid_argument("Can't initialize vector of dimensions != 3.");
        }
        _position = pos;
        _velocity = vel;
        _force = force;
        _mass = m;
        _kSpring = k;
        _cRepulsive = repulsive;
    }


    // set and get methods
    void setX(double x) {_position[0] = x;};
    void setY(double y) {_position[1] = y;};
    void setZ(double z) {_position[2] = z;};
    void setFx(double fx) {_force[0] = fx;};
    void setFy(double fy) {_force[1] = fy;};
    void setFz(double fz) {_force[2] = fz;};
    void addFx(double fx) {_force[0] += fx;};
    void addFy(double fy) {_force[1] += fy;};
    void addFz(double fz) {_force[2] += fz;};
    void setVx(double vx) {_velocity[0] = vx;};
    void setVy(double vy) {_velocity[1] = vy;};
    void setVz(double vz) {_velocity[2] = vz;};

    double getX() const {return _position[0];};
    double getY() const {return _position[1];};
    double getZ() const {return _position[2];};
    double getFx() const {return _force[0];};
    double getFy() const {return _force[1];};
    double getFz() const {return _force[2];};
    double getVx() const {return _velocity[0];};
    double getVy() const {return _velocity[1];};
    double getVz() const {return _velocity[2];};

    void setMass(double m) {_mass = m;};
    void setSpring(double k) {_kSpring = k;};
    void setRepulsive(double c) {_cRepulsive = c;};
    void setType(std::string name) {_type = name;};
    void setEquilibrium(double d) {_dEquilibrium = d;};

    double getMass() const {return _mass;};
    double getSpring() const {return _kSpring;};
    double getRepulsive() const {return _cRepulsive;};
    std::string getType() const {return _type;};
    double getEquilibrium() const {return _dEquilibrium;};

    void appendNeighbour(size_t i){
      _cutoff.push_back(i);
    };
    void appendSkin(size_t i){
      _skin.push_back(i);
    };
    int getCutoffN() {return _cutoff.size();};
    int getSkinN() {return _skin.size();};
    int getNeighbour(int i) {return _cutoff[i];};
    int getSkin(int i) {return _skin[i];};

    private:

      std::array<double,3> _position = {0.,0.,0.};
      std::array<double,3> _velocity = {0.,0.,0.};
      std::array<double,3> _force = {0.,0.,0.};

      double _mass, _kSpring, _cRepulsive, _dEquilibrium;
      std::string _type;

      // neighbours list. Append ids of the other atoms in the sphere
      std::vector<size_t> _cutoff; // vectors allow mutable size
      std::vector<size_t> _skin; // neighbours inside skin -> only forces ?

};

class MolecularDynamics : public Atom {

  public:

    // constructors
    MolecularDynamics(int n, std::vector<Atom> &atoms){
      _n = n;
      _mdAtoms = atoms;
    }

    // DISTANCE methods
    double distX(int i, int j){
      double delta = _mdAtoms[i].getX() - _mdAtoms[j].getX();
      if(_kbox > 0.) return delta;
      else return delta - round(delta/_side)*_side;
    };
    double distY(int i, int j){
      double delta = _mdAtoms[i].getY() - _mdAtoms[j].getY();
      if(_kbox > 0.) return delta;
      else return delta - round(delta/_side)*_side;
    };
    double distZ(int i, int j){
      double delta = _mdAtoms[i].getZ() - _mdAtoms[j].getZ();
      if(_kbox > 0.) return delta;
      else return delta - round(delta/_side)*_side;
    };
    double dist(int i, int j){
      return sqrt( square(distX(i,j)) + square(distY(i,j)) + square(distZ(i,j)) );
    };

    double distY(int i, int j){
      double delta = _mdAtoms[i].getY() - _mdAtoms[j].getY();
      if(_kbox > 0.) return delta;
      else return delta - round(delta/_side)*_side;
    };
    double distZ(int i, int j){
      double delta = _mdAtoms[i].getZ() - _mdAtoms[j].getZ();
      if(_kbox > 0.) return delta;
      else return delta - round(delta/_side)*_side;
    };
    double dist(int i, int j){
      return sqrt( square(distX(i,j)) + square(distY(i,j)) + square(distZ(i,j)) );
    };

    // NEIGHBOURS
    void calculateNeighbours(){
      double d = 0.;
      double limit = _cutoffLimit + _skindepth;
      // loop over all atoms
      for(size_t i=0; i<_n; i++){
        for(size_t j=i+1; j<_n; j++){
          d = dist(i,j);
          if(d>limit) continue;
          else if(d<_cutoffLimit){
            // append to inner neighbours
            _mdAtoms[i].appendNeighbour(j);
            _mdAtoms[j].appendNeighbour(i);
          } 
          else if(d<limit){
            // append to outer neighbours
            _mdAtoms[i].appendSkin(j);
            _mdAtoms[j].appendSkin(i);
          }
        }
      }

    }; 

    // FORCES
    void forceHelper(size_t i, size_t j){
      double f;
      double d=dist(i,j);
      double dx=distX(i,j);
      double dy=distY(i,j);
      double dz=distZ(i,j);
      if(_mdAtoms[i].getType() == _chainType){
        // chain atoms
        if(_mdAtoms[j].getType() == _chainType){
          f = - _mdAtoms[i].getSpring()*(d - _mdAtoms[i].getEquilibrium())/d;
          _mdAtoms[i].addFx(f*dx);
          _mdAtoms[i].addFy(f*dy);
          _mdAtoms[i].addFz(f*dz);
        }
      } else {
        // scattering
        f = _mdAtoms[i].getRepulsive()/pow(d,12);
        _mdAtoms[i].addFx(f*dx);
        _mdAtoms[i].addFy(f*dy);
        _mdAtoms[i].addFz(f*dz);
      }

      // box if no PBC
      if(_kbox>0.){
        double half = _side/2.;
        if(_mdAtoms[i].getX() < -half) _mdAtoms[i].addFx(-_kbox*(_mdAtoms[i].getX()+half));
        if(_mdAtoms[i].getX() >  half) _mdAtoms[i].addFx(-_kbox*(_mdAtoms[i].getX()-half));
        if(_mdAtoms[i].getY() < -half) _mdAtoms[i].addFy(-_kbox*(_mdAtoms[i].getY()+half));
        if(_mdAtoms[i].getY() >  half) _mdAtoms[i].addFy(-_kbox*(_mdAtoms[i].getY()-half));
        if(_mdAtoms[i].getZ() < -half) _mdAtoms[i].addFz(-_kbox*(_mdAtoms[i].getZ()+half));
        if(_mdAtoms[i].getZ() >  half) _mdAtoms[i].addFz(-_kbox*(_mdAtoms[i].getZ()-half));
      }

    };

    void calculateForces(){
      double f, d, dx,dy,dz;
      // loop over atoms INSIDE neighbours list
      for(size_t i=0; i<_n; i++){
        size_t nn = _mdAtoms[i].getCutoffN() + _mdAtoms[i].getSkinN();
        for(size_t m=0; m<nn; m++){
          int j = _mdAtoms[i].getNeighbour(m);
          forceHelper(i,j);          
        }
      }
    };

    // ENERGIES

    // dynamics integrators



  private:

    size_t _n; // number of total atoms
    std::vector<Atom> _mdAtoms;
    std::string _chainType;
    std::string _inputxyz, _paramfile, _trajfile;

    // parameters
    double _cutoffLimit, _skindepth;

    // boundary conditions stuff
    double _kbox, _side;
    int _dt, _nsteps;
    std::string _chaintype; // type of atoms in the molecular chain
    // TODO: this is allows for single atom-type chains, not more complex structures





}



#endif //ATOMS_H