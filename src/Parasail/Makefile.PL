use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Parasail',
    VERSION_FROM      => 'lib/Parasail.pm',
    LIBS              => ['-L. -lparasail'],
    INC               => '-I.',
);
