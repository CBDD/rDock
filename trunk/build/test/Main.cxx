#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include "RbtConfig.h"

const RbtString EXEVERSION = " ($Id$)";

int main( int argc, char **argv) {
  RbtString strExeName(argv[0]);
  Rbt::PrintStdHeader(std::cout,strExeName + EXEVERSION);
  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );
  runner.run();
  _RBTOBJECTCOUNTER_DUMP_(cout)
  return 0;
}
