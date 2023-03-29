package ESPRESSO_Version;

use strict;
use warnings;

use File::Basename qw(dirname);
use File::Temp qw(tempfile);
use IPC::Open3;

my $version_number = '1.3.1';

sub get_version_number {
    return $version_number;
}

sub get_source_code_commit {
    my $file_dir = dirname(__FILE__);
    my $temp_file_handle = tempfile();
    # Run command with stdout and stderr sent to $temp_file_handle
    my $git_pid = open3(my $unused_input_handle, $temp_file_handle,
                        $temp_file_handle,
                        "cd $file_dir; git log --pretty=format:%H -n 1");
    # If there was an error then report the commit as 'unknown'
    waitpid($git_pid, 0);
    if ($?) {
        close($temp_file_handle);
        return 'unknown';
    }
    # Read the commit info written by 'git log'
    my $git_output = readline($temp_file_handle);
    close($temp_file_handle);
    return $git_output;
}

sub parse_version_string {
  my $string = $_[0];
  my @parts = split /\./, $string;
  return @parts;
}

sub compare_versions {
  my ($parts_1_ref, $parts_2_ref) = @{$_[0]};
  my $len_1 = @{$parts_1_ref};
  my $len_2 = @{$parts_2_ref};
  my $shorter_len = $len_1 < $len_2 ? $len_1 : $len_2;
  for my $part_i (0 .. ($shorter_len - 1)) {
    my $part_1 = $parts_1_ref->[$part_i];
    my $part_2 = $parts_2_ref->[$part_i];
    if ($part_1 < $part_2) {
      return -1;
    }
    if ($part_2 < $part_1) {
      return 1;
    }
  }

  if ($len_1 < $len_2) {
    return -1
  }
  if ($len_2 < $len_1) {
    return 1;
  }
  return 0;
}

sub check_version_at_least {
  my ($parsed_version_ref, $target_string) = @{$_[0]};
  my @parsed_target = parse_version_string($target_string);
  if (@parsed_target == 0) {
    return 0;
  }
  my $compared = compare_versions([$parsed_version_ref, \@parsed_target]);
  return $compared >= 0;
}

sub check_samtools_version {
  my $version_command = 'samtools --version';
  my @samtools_version_lines = qx($version_command);
  if (@samtools_version_lines == 0) {
    return "No samtools version found from command: $version_command";
  }
  my @first_line_parts = split /\s+/, $samtools_version_lines[0];
  my $version_string = $first_line_parts[1];
  my @parsed_version = parse_version_string($version_string);
  if (@parsed_version == 0) {
    return "Could not parse version number from ($version_command): @samtools_version_lines";
  }
  my $required_version = '1.6';
  if (!check_version_at_least([\@parsed_version, $required_version])) {
    return "samtools must be at least version ($required_version), but found ($version_string)";
  }
  return "";
}

sub check_blast_version {
  my $version_command = 'blastn -version';
  my @blast_version_lines = qx($version_command);
  if (@blast_version_lines == 0) {
    return "No blastn version found from command: $version_command";
  }
  my @first_line_parts = split /\s+/, $blast_version_lines[0];
  my $version_string = $first_line_parts[1];
  my @parsed_version = parse_version_string($version_string);
  if (@parsed_version == 0) {
    return "Could not parse version number from ($version_command): @blast_version_lines";
  }
  my $required_version = '2.8.1';
  if (!check_version_at_least([\@parsed_version, $required_version])) {
    return "blastn must be at least version ($required_version), but found ($version_string)";
  }
  return "";
}

sub check_nhmmer_version {
  my $version_command = 'nhmmer -h';
  my @nhmmer_version_lines = qx($version_command);
  if (@nhmmer_version_lines < 2) {
    return "No nhmmer version found from command: $version_command";
  }
  my @second_line_parts = split /\s+/, $nhmmer_version_lines[1];
  my $version_string = $second_line_parts[2];
  my @parsed_version = parse_version_string($version_string);
  if (@parsed_version == 0) {
    return "Could not parse version number from ($version_command): @nhmmer_version_lines";
  }
  my $required_version = '3.3.1';
  if (!check_version_at_least([\@parsed_version, $required_version])) {
    return "nhmmer must be at least version ($required_version), but found ($version_string)";
  }
  return "";
}

# package must return true
1;
