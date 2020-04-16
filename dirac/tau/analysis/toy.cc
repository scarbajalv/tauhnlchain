#include <iostream>
#include <cmath>


using namespace std;

int main(int argc, char* argv[] ){

  float x = stof(argv[1]);

  switch (int(x-0.0001)){

    case -12: cout << "-12" << endl; break;
    case -14: cout << "-14" << endl; break;
    case -16: cout << "-16" << endl; break;
    default: break;

  }

  return 0;

}