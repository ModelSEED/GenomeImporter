package GenomeImporter::GenomeImporterImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = '0.0.1';
our $GIT_URL = '';
our $GIT_COMMIT_HASH = '';

=head1 NAME

GenomeImporter

=head1 DESCRIPTION

A KBase module: GenomeImporter
This sample module contains one small method - filter_contigs.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use AssemblyUtil::AssemblyUtilClient;
use KBaseReport::KBaseReportClient;
use Bio::ModelSEED::Client::SAP;
use Bio::KBase::workspace::Client;
use Bio::KBase::utilities;
use Bio::KBase::ObjectAPI::utilities;
use Bio::KBase::kbaseenv;
use GenomeAnnotationAPI::GenomeAnnotationAPIClient;
use AssemblyUtil::AssemblyUtilClient;
use Data::UUID;

sub get_PATRIC_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["id","source","workspace"],{});
	my $genomeid = $args->{id};
	my $source = $args->{source};
	my $refseq = 0;
	if ($source eq "patricrefseq") {
		$refseq = 1;
	}
	my $data = Bio::KBase::ObjectAPI::utilities::rest_download({url => Bio::KBase::utilities::conf("GenomeImporter","data_api_url")."genome/?genome_id=".$genomeid."&http_accept=application/json",token => Bio::KBase::utilities::token()});
	$data = $data->[0];
	$data = Bio::KBase::utilities::args($data,[],{
		genome_length => 0,
		contigs => 0,
		genome_name => "Unknown",
		taxon_lineage_names => ["Unknown"],
		owner => "Unknown",
		gc_content => 0,
		publication => "Unknown",
		completion_date => "1970-01-01T00:00:00+0000"
	});
    my $genomesource = "PATRIC";
    if ($refseq == 1) {
    	$genomesource = "PATRICRefSeq";
    }
	my $genome = {
    	id => $genomeid,
		scientific_name => $data->{genome_name},
		domain => $data->{taxon_lineage_names}->[0],
		genetic_code => 11,
		dna_size => $data->{genome_length},
		num_contigs => 0,
		contig_lengths => [],
		contig_ids => [],
		source => $genomesource,
		source_id => $genomeid,
		md5 => "none",
		taxonomy => join(":",@{$data->{taxon_lineage_names}}),
		gc_content => $data->{gc_content},
		complete => 1,
		publications => [$data->{publication}],
		features => []
	};
	#Retrieving feature information
	my $start = 0;
	my $params = {};
	my $loopcount = 0;
	my $ftrcount = 0;
	my $allftrs = [];
	while ($start >= 0 && $loopcount < 100) {
		$loopcount++;#Insurance that no matter what, this loop won't run for more than 100 iterations
		my $ftrdata = Bio::KBase::ObjectAPI::utilities::rest_download({url => Bio::KBase::utilities::conf("GenomeImporter","data_api_url")."genome_feature/?genome_id=".$genomeid."&http_accept=application/json&limit(10000,$start)",token => Bio::KBase::utilities::token()},$params);
		if (defined($ftrdata) && @{$ftrdata} > 0) {
			push(@{$allftrs},@{$ftrdata});
		}
		my $currentcount = @{$ftrdata};
		$ftrcount += $currentcount;
		if ($ftrcount < $params->{count}) {
			$start = $ftrcount;
		} else {
			$start = -1;
		}
	}
	my $patricids = {};
	my $refseqids = {};
	my $stops = {};
	for (my $i=0; $i < @{$allftrs}; $i++) {
		my $ftrdata = $allftrs->[$i];
		if ($ftrdata->{feature_type} ne "pseudogene") {
			if ($ftrdata->{strand} eq "-" && $ftrdata->{annotation} eq "RefSeq") {
				$stops->{$ftrdata->{start}} = $ftrdata;
			} elsif ($ftrdata->{annotation} eq "RefSeq") {
				$stops->{$ftrdata->{end}} = $ftrdata;
			}
		}
	}
	my $refseqgenes = 0;
	my $patricgenes = 0;
	my $match = 0;
	my $weakermatch = 0;
	my $unsortedftrlist = [];
	for (my $i=0; $i < @{$allftrs}; $i++) {
		my $ftrdata = $allftrs->[$i];
		if ($ftrdata->{annotation} eq "PATRIC") {
			push(@{$unsortedftrlist},$ftrdata);
			$patricgenes++;
			if ($ftrdata->{strand} eq "-" && defined($stops->{$ftrdata->{start}}) && $ftrdata->{strand} eq $stops->{$ftrdata->{start}}->{strand}){
				$match++;
				$ftrdata->{refseqgene} = $stops->{$ftrdata->{start}};
				if (defined($stops->{$ftrdata->{start}}->{patricgene})) {
					delete($stops->{$ftrdata->{start}}->{refseqgene});
					$weakermatch--;
				}
				$stops->{$ftrdata->{start}}->{patricgene} = $ftrdata->{feature_id};
			} elsif ($ftrdata->{strand} eq "+" && defined($stops->{$ftrdata->{end}}) && $ftrdata->{strand} eq $stops->{$ftrdata->{end}}->{strand}){
				$match++;
				$ftrdata->{refseqgene} = $stops->{$ftrdata->{end}};
				if (defined($stops->{$ftrdata->{end}}->{patricgene})) {
					delete($stops->{$ftrdata->{end}}->{refseqgene});
					$weakermatch--;
				}
				$stops->{$ftrdata->{end}}->{patricgene} = $ftrdata->{feature_id};
#			} elsif ($ftrdata->{strand} eq "+") {
#				for (my $j=0; $j < 100; $j++) {
#					my $startindex = ($ftrdata->{end} - 50 + $j);
#					if (defined($stops->{$startindex}) && !defined($stops->{$startindex}->{patricgene}) && $ftrdata->{strand} eq $stops->{$startindex}->{strand} && abs($ftrdata->{start}-$stops->{$startindex}->{start}) < 50) {
#						$weakermatch++;
#						$ftrdata->{refseqgene} = $stops->{$startindex};
#						$stops->{$startindex}->{patricgene} = $ftrdata->{feature_id};
#						last;
#					}
#				}
#			} elsif ($ftrdata->{strand} eq "-") {
#				for (my $j=0; $j < 100; $j++) {
#					my $startindex = ($ftrdata->{start} - 50 + $j);
#					if (defined($stops->{$startindex}) && !defined($stops->{$startindex}->{patricgene}) && $ftrdata->{strand} eq $stops->{$startindex}->{strand} && abs($ftrdata->{end}-$stops->{$startindex}->{end}) < 50) {
#						$weakermatch++;
#						$ftrdata->{refseqgene} = $stops->{$startindex};
#						$stops->{$startindex}->{patricgene} = $ftrdata->{feature_id};
#						last;
#					}
#				}
			}
		} else {
			$refseqgenes++;
			$refseqids->{$ftrdata->{feature_id}} = $ftrdata;
		}
	}
	my $sortedftrlist = [sort { $b->{start} cmp $a->{start} } @{$unsortedftrlist}];
	my $funchash = Bio::KBase::ObjectAPI::utilities::get_SSO();
	for (my $i=0; $i < @{$sortedftrlist}; $i++) {
		my $id;
		my $data = $sortedftrlist->[$i];
		if ($refseq == 1 && defined($data->{refseqgene}->{refseq_locus_tag})) {
			$id = $data->{refseqgene}->{refseq_locus_tag}; 
		} else {
			$id = $data->{patric_id};
		}
		if (defined($id)) {		
			my $ftrobj = {id => $id,type => "CDS",aliases=>[]};		
			if (defined($data->{start})) {
				$ftrobj->{location} = [[$data->{sequence_id},$data->{start},$data->{strand},$data->{na_length}]];
			}
			if (defined($data->{feature_type})) {
				$ftrobj->{type} = $data->{feature_type};
			}
			if (defined($data->{product})) {
				$ftrobj->{function} = $data->{product};
			}
			if (defined($data->{na_sequence})) {
				$ftrobj->{dna_sequence} = $data->{na_sequence};
				$ftrobj->{dna_sequence_length} = $data->{na_length};
			}
			if (defined($data->{aa_sequence})) {
				$ftrobj->{protein_translation} = $data->{aa_sequence};
				$ftrobj->{protein_translation_length} = $data->{aa_length};
				$ftrobj->{md5} = $data->{aa_sequence_md5};
			}
			if (defined($data->{refseqgene}->{refseq_locus_tag})) {
				$ftrobj->{ontology_terms}->{RefSeq}->{$data->{refseqgene}->{refseq_locus_tag}} = {
					id => $data->{refseqgene}->{refseq_locus_tag},
					term_name => $data->{refseqgene}->{product}
				};
				if (defined($data->{refseqgene}->{gene})) {
					$ftrobj->{ontology_terms}->{GeneName}->{$data->{refseqgene}->{gene}} = {
						id => $data->{refseqgene}->{gene}
					};
				}
			}
			if (defined($data->{patric_id})) {
				$ftrobj->{ontology_terms}->{PATRIC}->{$data->{patric_id}} = {
					id => $data->{patric_id}
				};
			}
			if (defined($data->{figfam_id})) {
				$ftrobj->{ontology_terms}->{FigFam}->{$data->{figfam_id}} = {
					id => $data->{figfam_id}
				};
			}
			if (defined($data->{pgfam_id})) {
				$ftrobj->{ontology_terms}->{PGFam}->{$data->{pgfam_id}} = {
					id => $data->{pgfam_id}
				};
			}
			if (defined($data->{plfam_id})) {
				$ftrobj->{ontology_terms}->{PLFam}->{$data->{plfam_id}} = {
					id => $data->{plfam_id}
				};
			}
			if (defined($data->{accession})) {
				$ftrobj->{ontology_terms}->{Accession}->{$data->{accession}} = {
					id => $data->{accession}
				};
			}
			if (defined($data->{product})) {
				my $function = $data->{product};
	  			my $array = [split(/\#/,$function)];
	  			$function = shift(@{$array});
				$function =~ s/\s+$//;
				$array = [split(/\s*;\s+|\s+[\@\/]\s+/,$function)];
				for (my $k=0; $k < @{$array}; $k++) {
					my $rolename = lc($array->[$k]);
					$rolename =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
					$rolename =~ s/\s//g;
					$rolename =~ s/\#.*$//g;
					if (defined($funchash->{$rolename})) {
						$ftrobj->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}} = {
							id => $funchash->{$rolename}->{id},
							term_name => $funchash->{$rolename}->{name}
						};
					}
				}
			}
			if (defined($data->{go})) {
				for (my $k=0; $k < @{$data->{go}}; $k++) {
					my $array = [split(/\|/,$data->{go}->[$k])];
					$ftrobj->{ontology_terms}->{GO}->{$array->[0]} = {
						id => $array->[0],
						term_name => $array->[1]
					};
				}
			}
			push(@{$genome->{features}},$ftrobj);
		}
	}
    $data = Bio::KBase::ObjectAPI::utilities::rest_download({url => Bio::KBase::utilities::conf("GenomeImporter","data_api_url")."genome_sequence/?genome_id=".$genomeid."&http_accept=application/json",token => Bio::KBase::utilities::token()});
	my $contigHash;
	my $contigarray;
	for (my $i=0; $i < @{$data}; $i++) {
		$contigHash->{$data->[$i]->{sequence_id}} = $data->[$i]->{sequence};
		push(@{$contigarray},$data->[$i]->{sequence});
		$genome->{num_contigs}++;
		push(@{$genome->{contig_lengths}},length($data->[$i]->{sequence}));
		push(@{$genome->{contig_ids}},$data->[$i]->{sequence_id});
	}
	my $sortedcontigs = [sort { $a cmp $b } @{$contigarray}];
	my $str = "";
	for (my $i=0; $i < @{$sortedcontigs}; $i++) {
		if (length($str) > 0) {
			$str .= ";";
		}
		$str .= $sortedcontigs->[$i];	
	}
	$genome->{md5} = Digest::MD5::md5_hex($str);
	$genome->{assembly_ref} = $self->save_assembly({
		workspace => $args->{workspace},
		data => $contigHash,
		name => $genomeid.".contigs"
	});
	$self->save_genome({
		workspace => $args->{workspace},
		data => $genome,
		name => $genomeid
	});
}

sub get_SEED_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["id","source","workspace"],{});
	my $id = $args->{id};
	my $source = $args->{source};
	my $sapsvr;
	if ($source eq "pubseed") {
		$sapsvr = Bio::ModelSEED::Client::SAP->new();
	} elsif ($source eq "coreseed") {
		$sapsvr = Bio::ModelSEED::Client::SAP->new({url => "http://core.theseed.org/FIG/sap_server.cgi"});
	}
	my $translation = {
		pubseed => "PubSEED",
		coreseed => "CoreSEED"
	};
	my $data = $sapsvr->genome_data({
		-ids => [$id],
		-data => [qw(gc-content dna-size name taxonomy domain genetic-code)]
	});
	if (!defined($data->{$id})) {
    	$self->error("PubSEED genome ".$id." not found!");
    }
    my $genomeObj = {
		id => $id,
		scientific_name => $data->{$id}->[2],
		domain => $data->{$id}->[4],
		genetic_code => $data->{$id}->[5],
		dna_size => $data->{$id}->[1],
		num_contigs => 0,
		contig_lengths => [],
		contig_ids => [],
		source => $translation->{$source},
		source_id => $id,
		taxonomy => $data->{$id}->[3],
		gc_content => $data->{$id}->[0]/100,
		complete => 1,
		publications => [],
		features => [],
    };
	my $featureHash = $sapsvr->all_features({-ids => $id});
	my $genomeHash = $sapsvr->genome_contigs({
		-ids => [$id]
	});
	my $featureList = $featureHash->{$id};
	my $contigList = $genomeHash->{$id};
	my $functions = $sapsvr->ids_to_functions({-ids => $featureList});
	my $dnaseq = $sapsvr->ids_to_sequences({-ids => $featureList});
	
	my $locations = $sapsvr->fid_locations({-ids => $featureList});
	my $sequences = $sapsvr->fids_to_proteins({-ids => $featureList,-sequence => 1});
	my $contigHash = $sapsvr->contig_sequences({
		-ids => $contigList
	});
	my $contigarray;
	foreach my $contigid (keys(%{$contigHash})) {
		push(@{$contigarray},$contigHash->{$contigid});
		$genomeObj->{num_contigs}++;
		push(@{$genomeObj->{contig_lengths}},length($contigHash->{$contigid}));
		push(@{$genomeObj->{contig_ids}},$contigid);
	}
	my $sortedcontigs = [sort { $a cmp $b } @{$contigarray}];
	my $str = "";
	for (my $i=0; $i < @{$sortedcontigs}; $i++) {
		if (length($str) > 0) {
			$str .= ";";
		}
		$str .= $sortedcontigs->[$i];	
	}
	$genomeObj->{md5} = Digest::MD5::md5_hex($str);
	for (my $i=0; $i < @{$featureList}; $i++) {
		my $feature = {
  			id => $featureList->[$i],
			type => "peg",
			publications => [],
			subsystems => [],
			protein_families => [],
			aliases => [],
			annotations => [],
			subsystem_data => [],
			regulon_data => [],
			atomic_regulons => [],
			coexpressed_fids => [],
			co_occurring_fids => []
  		};
  		if ($featureList->[$i] =~ m/\.([^\.]+)\.\d+$/) {
  			$feature->{type} = $1;
  		}
		if (defined($functions->{$featureList->[$i]})) {
			$feature->{function} = $functions->{$featureList->[$i]};
		}
		if (defined($sequences->{$featureList->[$i]})) {
			$feature->{protein_translation} = $sequences->{$featureList->[$i]};
			$feature->{protein_translation_length} = length($feature->{protein_translation});
  			$feature->{dna_sequence_length} = 3*$feature->{protein_translation_length};
  			$feature->{md5} = Digest::MD5::md5_hex($feature->{protein_translation});
		}
		if (defined($dnaseq->{$featureList->[$i]})) {
			$feature->{dna_sequence} = $dnaseq->{$featureList->[$i]};
			$feature->{dna_sequence_length} = length($dnaseq->{$featureList->[$i]});
		}
  		if (defined($locations->{$featureList->[$i]}->[0])) {
			for (my $j=0; $j < @{$locations->{$featureList->[$i]}}; $j++) {
				my $loc = $locations->{$featureList->[$i]}->[$j];
				if ($loc =~ m/^(.+)_(\d+)([\+\-])(\d+)$/) {
					my $array = [split(/:/,$1)];
					if ($3 eq "-" || $3 eq "+") {
						$feature->{location}->[$j] = [$array->[1],$2,$3,$4];
					} elsif ($2 > $4) {
						$feature->{location}->[$j] = [$array->[1],$2,"-",($2-$4)];
					} else {
						$feature->{location}->[$j] = [$array->[1],$2,"+",($4-$2)];
					}
					$feature->{location}->[$j]->[1] = $feature->{location}->[$j]->[1]+0;
					$feature->{location}->[$j]->[3] = $feature->{location}->[$j]->[3]+0;
				}
			}
			
		}
  		push(@{$genomeObj->{features}},$feature);	
	}
	$genomeObj->{assembly_ref} = $self->save_assembly({
		workspace => $args->{workspace},
		data => $contigHash,
		name => $id.".contigs"
	});
	$self->save_genome({
		workspace => $args->{workspace},
		data => $genomeObj,
		name => $id
	});
}

sub save_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["workspace","data"],{
		name => $args->{data}->{id}
	});
	my $ga = new GenomeAnnotationAPI::GenomeAnnotationAPIClient(Bio::KBase::utilities::conf("GenomeImporter","call_back_url"));
	my $gaout = $ga->save_one_genome_v1({
		workspace => $args->{workspace},
        name => $args->{name},
        data => $args->{data},
        provenance => [],
        hidden => 0
	});
	Bio::KBase::kbaseenv::add_object_created({
		"ref" => $gaout->{info}->[6]."/".$gaout->{info}->[0]."/".$gaout->{info}->[4],
		"description" => "Genome object for ".$args->{name}
	});
	return $gaout->{info}->[6]."/".$gaout->{info}->[0]."/".$gaout->{info}->[4];
}

sub save_assembly {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["workspace","data","name"],{});
	my $filename = Bio::KBase::utilities::conf("GenomeImporter","scratch")."/".$args->{name}.".fa";
	open(my $fo, ">", $filename) || die "Could not open file ".$filename;
	foreach my $contigid (keys(%{$args->{data}})) {
		print $fo "<".$contigid."\n";
		my $sequence = $args->{data}->{$contigid};
		my $seq_len = length($sequence);
		for (my $i = 0; $i < $seq_len; $i += 60) {
		    my $segment = substr($sequence, $i, 60);
		    print $fo $segment."\n";
		}
	}
	close($fo);
	my $ga = new AssemblyUtil::AssemblyUtilClient(Bio::KBase::utilities::conf("GenomeImporter","call_back_url"));
	my $assemblyref = $ga->save_assembly_from_fasta({
		file => $filename,
        workspace_name => $args->{workspace},
        assembly_name => $args->{name}
	});
	Bio::KBase::kbaseenv::add_object_created({
		"ref" => $assemblyref,
		"description" => "Contigs for ".$args->{name}
	});
	return $assemblyref;
}

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    Bio::KBase::utilities::read_config({
		filename => $ENV{KB_DEPLOYMENT_CONFIG},
		service => "GenomeImporter"
	});
    Bio::KBase::utilities::setconf("GenomeImporter","call_back_url",$ENV{ SDK_CALLBACK_URL });
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 import_external_genome

  $output = $obj->import_external_genome($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a GenomeImporter.ImportGenomeParams
$output is a GenomeImporter.ImportGenomeResults
ImportGenomeParams is a reference to a hash where the following keys are defined:
	genome_ids has a value which is a string
	source has a value which is a string
	workspace has a value which is a string
ImportGenomeResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a GenomeImporter.ImportGenomeParams
$output is a GenomeImporter.ImportGenomeResults
ImportGenomeParams is a reference to a hash where the following keys are defined:
	genome_ids has a value which is a string
	source has a value which is a string
	workspace has a value which is a string
ImportGenomeResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text



=item Description

Function to import a list of genomes from a specified source

=back

=cut

sub import_external_genome
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to import_external_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'import_external_genome');
    }

    my $ctx = $GenomeImporter::GenomeImporterServer::CallContext;
    my($output);
    #BEGIN import_external_genome
    my($args) = @_;
	$args = Bio::KBase::utilities::args($args,["genome_ids","workspace","source"],{});
    my $genomes = [split(/[\n;\|]+/,$args->{genome_ids})];
    for (my $i=0; $i<@{$genomes};$i++) {
    	print "Now importing ".$genomes->[$i]." from ".$args->{source}."\n";
    	if ($args->{source} eq "pubseed" || $args->{source} eq "coreseed") {
    		my $refs = $self->get_SEED_genome({
    			id => $genomes->[$i],
    			source => $args->{source},
    			workspace => $args->{workspace}
    		});
    	} elsif ($args->{source} eq "patric" || $args->{source} eq "patricrefseq") {
    		my $refs = $self->get_PATRIC_genome({
    			id => $genomes->[$i],
    			source => $args->{source},
    			workspace => $args->{workspace}
    		});
    	}
    }
    my $reportout = Bio::KBase::kbaseenv::create_report({
    	workspace_name => $args->{workspace},
    	report_object_name => Bio::KBase::utilities::processid()
    });
    $output->{report_ref} = $reportout->{"ref"};
	$output->{report_name} = Bio::KBase::utilities::processid();
    #END import_external_genome
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to import_external_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'import_external_genome');
    }
    return($output);
}




=head2 status 

  $return = $obj->status()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module status. This is a structure including Semantic Versioning number, state and git info.

=back

=cut

sub status {
    my($return);
    #BEGIN_STATUS
    $return = {"state" => "OK", "message" => "", "version" => $VERSION,
               "git_url" => $GIT_URL, "git_commit_hash" => $GIT_COMMIT_HASH};
    #END_STATUS
    return($return);
}

=head1 TYPES



=head2 ImportGenomeParams

=over 4



=item Description

Input parameters for the import_external_genome function


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_ids has a value which is a string
source has a value which is a string
workspace has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_ids has a value which is a string
source has a value which is a string
workspace has a value which is a string


=end text

=back



=head2 ImportGenomeResults

=over 4



=item Description

Output structure for the import_external_genome function


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=cut

1;
