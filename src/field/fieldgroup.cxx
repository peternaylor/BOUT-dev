
#include <bout/fieldgroup.hxx>

FieldGroup::FieldGroup(FieldData &f) {
  fvec.push_back(&f);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2) {
  fvec.push_back(&f1); fvec.push_back(&f2);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4); fvec.push_back(&f5);
}

FieldGroup::FieldGroup(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4, FieldData &f5, FieldData &f6) {
  fvec.push_back(&f1); fvec.push_back(&f2); fvec.push_back(&f3); fvec.push_back(&f4);
  fvec.push_back(&f5); fvec.push_back(&f6);
}

FieldGroup FieldGroup::operator+(const FieldGroup &other){
  FieldGroup temp=(*this); //Temporary field group -- Initialise temporary to hold contents of this

  //Now add contents of other
  for(int i=0;i<other.fvec.size();i++){
    temp.add(*(other.fvec[i]));
  };

  //Return copy of temp
  return temp;
};

FieldGroup& FieldGroup::operator+=(const FieldGroup &other){
  (*this)=(*this)+other;
  return *this;
};
