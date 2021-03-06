/**************************************************************************
 * Generate a field with specified values, mainly for creating
 * initial perturbations
 * 
 * 
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/


class FieldFactory;

#ifndef __FIELD_FACTORY_H__
#define __FIELD_FACTORY_H__

#include "bout/mesh.hxx"

#include "bout/sys/expressionparser.hxx"

#include "field2d.hxx"
#include "field3d.hxx"
#include "options.hxx"

#include <string>
#include <map>
#include <list>

// Utility routines to create generators from values

FieldGenerator* generator(BoutReal value);
FieldGenerator* generator(BoutReal *ptr);

//////////////////////////////////////////////////////////
// Create a tree of generators from an input string

class FieldFactory : public ExpressionParser {
public:
  FieldFactory(Mesh *m, Options *opt = NULL);
  ~FieldFactory();
  
  const Field2D create2D(const std::string &value, Options *opt = NULL, Mesh *m = NULL, CELL_LOC loc=CELL_CENTRE, BoutReal t=0.0);
  const Field3D create3D(const std::string &value, Options *opt = NULL, Mesh *m = NULL, CELL_LOC loc=CELL_CENTRE, BoutReal t=0.0);

  // Parse a string into a tree of generators
  FieldGenerator* parse(const std::string &input, Options *opt=NULL);

  // Singleton object
  static FieldFactory *get();
protected:
  // These functions called by the parser
  FieldGenerator* resolve(std::string &name);
  
private:
  Mesh *fieldmesh;  
  Options *options;

  std::list<std::string> lookup; // Names currently being parsed
  
  // Cache parsed strings
  std::map<std::string, FieldGenerator* > cache;
  
  Options* findOption(Options *opt, const std::string &name, std::string &val);
};

//////////////////////////////////////////////////////////
// Generator from function

class FieldFunction : public FieldGenerator {
public:
  FieldFunction(FuncPtr userfunc) : func(userfunc) {}
  double generate(double x, double y, double z, double t) {
    return func(t, x, y, z);
  }
private:
  FieldFunction();
  
  FuncPtr func;
};

//////////////////////////////////////////////////////////
// Null generator 

class FieldNull : public FieldGenerator {
public:
  double generate(double x, double y, double z, double t) {
    return 0.0;
  }
  FieldGenerator* clone(const std::list<FieldGenerator*> args) {
    return this;
  }
  /// Singeton
  static FieldGenerator* get() {
    static FieldNull *instance = 0;
    
    if(!instance)
      instance = new FieldNull();
    return instance;
  }
private:
};

#endif // __FIELD_FACTORY_H__
