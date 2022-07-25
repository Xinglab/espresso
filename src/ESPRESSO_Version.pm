package ESPRESSO_Version;

use strict;
use warnings;

use File::Basename qw(dirname);

my $version_number = 'beta1.3.0';

sub get_version_number {
    return $version_number;
}

sub get_source_code_commit {
    my $file_dir = dirname(__FILE__);
    my $git_output = qx(cd $file_dir; git log --pretty=format:%H -n 1);
    if ($?) {
        return 'unknown';
    }
    return $git_output;
}

# package must return true
1;
