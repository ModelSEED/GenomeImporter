use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::utilities;
use Bio::KBase::kbaseenv;
use GenomeImporter::GenomeImporterImpl;

#my $tester = LocalTester->new($ENV{'KB_DEPLOYMENT_CONFIG'});
my $tester = LocalTester->new("/Users/chenry/code/GenomeImporter/localdeploy.cfg");
$tester->run_tests();

{
	package LocalTester;
	use strict;
	use Test::More;
    sub new {
        my ($class,$configfile) = @_;
        Bio::KBase::kbaseenv::create_context_from_client_config({
        	filename => "/Users/chenry/.kbase_config"
        });
        my $c = Bio::KBase::utilities::read_config({
        	filename => $configfile,
			service => 'GenomeImporter'
        });
        my $object = GenomeImporter::GenomeImporterImpl->new();
        my $self = {
            token => Bio::KBase::utilities::token(),
            config_file => $configfile,
            config => $c->{fba_tools},
            user_id => Bio::KBase::utilities::user_id(),
            ws_client => Bio::KBase::kbaseenv::ws_client(),
            obj => $object,
            testcount => 0,
            completetestcount => 0,
            dumpoutput => 0,
            testoutput => {},
            showerrors => 1
        };
        return bless $self, $class;
    }
    sub test_harness {
		my($self,$function,$parameters,$name,$tests,$fail_to_pass,$dependency) = @_;
		$self->{testoutput}->{$name} = {
			output => undef,
			"index" => $self->{testcount},
			tests => $tests,
			command => $function,
			parameters => $parameters,
			dependency => $dependency,
			fail_to_pass => $fail_to_pass,
			pass => 1,
			function => 1,
			status => "Failed initial function test!"
		};
		$self->{testcount}++;
		if (defined($dependency) && $self->{testoutput}->{$dependency}->{function} != 1) {
			$self->{testoutput}->{$name}->{pass} = -1;
			$self->{testoutput}->{$name}->{function} = -1;
			$self->{testoutput}->{$name}->{status} = "Test skipped due to failed dependency!";
			return;
		}
		my $output;
		#eval {
			if (defined($parameters)) {
				$output = $self->{obj}->$function($parameters);
			} else {
				$output = $self->{obj}->$function();
			}
		#};
		my $errors;
		if ($@) {
			$errors = $@;
		}
		$self->{completetestcount}++;
		if (defined($output)) {
			$self->{testoutput}->{$name}->{output} = $output;
			$self->{testoutput}->{$name}->{function} = 1;
			if (defined($fail_to_pass) && $fail_to_pass == 1) {
				$self->{testoutput}->{$name}->{pass} = 0;
				$self->{testoutput}->{$name}->{status} = $name." worked, but should have failed!"; 
				ok $self->{testoutput}->{$name}->{pass} == 1, $self->{testoutput}->{$name}->{status};
			} else {
				ok 1, $name." worked as expected!";
				for (my $i=0; $i < @{$tests}; $i++) {
					$self->{completetestcount}++;
					$tests->[$i]->[2] = eval $tests->[$i]->[0];
					if ($tests->[$i]->[2] == 0) {
						$self->{testoutput}->{$name}->{pass} = 0;
						$self->{testoutput}->{$name}->{status} = $name." worked, but sub-tests failed!"; 
					}
					ok $tests->[$i]->[2] == 1, $tests->[$i]->[1];
				}
			}
		} else {
			$self->{testoutput}->{$name}->{function} = 0;
			if (defined($fail_to_pass) && $fail_to_pass == 1) {
				$self->{testoutput}->{$name}->{pass} = 1;
				$self->{testoutput}->{$name}->{status} = $name." failed as expected!";
			} else {
				$self->{testoutput}->{$name}->{pass} = 0;
				$self->{testoutput}->{$name}->{status} = $name." failed to function at all!";
			}
			ok $self->{testoutput}->{$name}->{pass} == 1, $self->{testoutput}->{$name}->{status};
			if ($self->{showerrors} && $self->{testoutput}->{$name}->{pass} == 0 && defined($errors)) {
				print "Errors:\n".$errors."\n";
			}
		}
		if ($self->{dumpoutput}) {
			print "$function output:\n".Data::Dumper->Dump([$output])."\n\n";
		}
		return $output;
	}
	sub run_tests {
		my($self) = @_;
		print "TEST1\n";
		my $output = $self->test_harness("import_external_genome",{
			genome_ids => "83333.1",
			workspace => "chenry:1491801790774",
			source => "pubseed"
		},"import pubseed genome",[],0,undef);
		print "TEST2\n";
		$output = $self->test_harness("import_external_genome",{
			genome_ids => "100226.1",
			workspace => "chenry:1491801790774",
			source => "coreseed"
		},"import pubseed genome",[],0,undef);
		print "TEST3\n";
		$output = $self->test_harness("import_external_genome",{
			genome_ids => "1110693.3",
			workspace => "chenry:1491801790774",
			source => "patric"
		},"import patric genome",[],0,undef);
		print "TEST4\n";
		$output = $self->test_harness("import_external_genome",{
			genome_ids => "1110693.3",
			workspace => "chenry:1491801790774",
			source => "patricrefseq"
		},"import patric genome",[],0,undef);
	}
}