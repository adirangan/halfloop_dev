system('make -f halfloop_OptiPlex.make clean; make -f halfloop_OptiPlex.make halfloop_lib;');
if (libisloaded('halfloop_lib')); unloadlibrary('halfloop_lib'); end;
loadlibrary('halfloop_lib','halfloop_lib.h','includepath','/home/rangan/dir_bcc/dir_halfloop_dev/');
libfunctions('halfloop_lib');
