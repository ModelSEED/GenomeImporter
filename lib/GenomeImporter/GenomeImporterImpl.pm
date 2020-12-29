package GenomeImporter::GenomeImporterImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = '0.0.1';
our $GIT_URL = 'ssh://git@github.com/ModelSEED/GenomeImporter';
our $GIT_COMMIT_HASH = '3d6d3dcb0e4d07a851b5e540154e0058e6094908';

=head1 NAME

GenomeImporter

=head1 DESCRIPTION

A KBase module: GenomeImporter
This sample module contains one small method - filter_contigs.

=cut

#BEGIN_HEADER
use Digest::MD5;
use KBaseReport::KBaseReportClient;
use Bio::ModelSEED::Client::SAP;
use Bio::KBase::utilities;
use Bio::KBase::ObjectAPI::utilities;
use Bio::KBase::kbaseenv;
use GenomeAnnotationAPI::GenomeAnnotationAPIClient;
use annotation_ontology_api::annotation_ontology_apiServiceClient;
use AssemblyUtil::AssemblyUtilClient;
use Data::UUID;
use P3DataAPI;

my $d = undef;

sub query_for_sequences {
	my($self,$md5s,$md5_hash,$type) = @_;
	if (!defined($d)) {
		$d = P3DataAPI->new();
	}
	my $batchsize = 100;
    my $n = @{$md5s};
    print $type." count:".$n."\n";
    my $end = 99;
    for (my $i = 0; $i < $n; $i = $end + 1) {
        $end = ($i + $batchsize) > $n ? ($n - 1) : ($i + $batchsize - 1);
        my $trycount = 0;
       	my $query = "(";
       	for (my $j=$i;$j < $end;$j++) {
       		if (length($query) > 2) {
       			$query .= ",";
       		}
       		$query .= $md5s->[$j];
       	}
       	$query .= ")";
        my $seqres;
        while (1) {
        		eval {
	        		$seqres = [$d->query('feature_sequence',
	                        ['select', 'sequence,md5,sequence_type'],
	                        ['in', 'md5', $query])]; 
	        };
	        if ($@ && $trycount < 3) {
	        		sleep(5);
				$trycount++;
			} elsif ($@) {
				$self->error("Too many failures trying to retrieve sequence data!");
				print $@;
			} else {
				last;
			}
        } 
		for my $seqent (@{$seqres}) {
	    		if (defined($md5_hash->{$seqent->{md5}})) {
	    			$md5_hash->{$seqent->{md5}}->{$type."_sequence"} = $seqent->{sequence};
	    			$md5_hash->{$seqent->{md5}}->{$type."_length"} = length($seqent->{sequence});
	    		} else {
	    			print $i." ".$end." MD5 not found:".$seqent->{md5}."\n";
	    		}
	    	}
    }
}

sub get_base_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["id","source","gc_content","assembly_ref","dna_size","md5","scientific_name","taxonomy"],{
		genetic_code => 11,
		contig_ids => [],
		contig_lengths => [],
		cdss => [],
		features => [],
		non_coding_features => [],
		mrnas => [],
		num_contigs => 0,
		external_source_origination_date => Bio::KBase::utilities::timestamp(1),
		notes => "Genome imported from ".$args->{source},
		source_id => $args->{id},
		domain => "Bacteria",
		taxon_ref => "ReferenceTaxons/unknown_taxon",
		feature_counts => {
	        CDS => 0,
	        gene => 0,
	        ncRNA => 0,
	        "non-protein_encoding_gene" => 0,
	        protein_encoding_gene => 0,
	        rRNA => 0,
	        regulatory => 0,
	        repeat_region => 0,
	        tRNA => 0,
	        tmRNA => 0
    		}
	});
	return {
    		id => $args->{id},
    		molecule_type => "DNA",
    		md5 => $args->{md5},
    		gc_content => $args->{gc_content}+0,
    		dna_size => $args->{dna_size}+0,
    		genetic_code => $args->{genetic_code}+0,
    		domain => "Bacteria",
    		genome_tiers => [
			"ExternalDB",
			"User"
		],
		publications => [],
		scientific_name => $args->{scientific_name},
		source => $args->{source},
		source_id => $args->{source_id},
		assembly_ref => $args->{assembly_ref},
		external_source_origination_date => $args->{external_source_origination_date},
		notes => $args->{notes},
    		taxon_ref => $args->{taxon_ref},
    		taxonomy => $args->{taxonomy},
    		num_contigs => $args->{num_contigs}+0,
		contig_ids => $args->{contig_ids},
		contig_lengths => $args->{contig_lengths},
		features => $args->{features},
		non_coding_features => $args->{non_coding_features},
		cdss => $args->{cdss},
		mrnas => $args->{mrnas},
		feature_counts => $args->{feature_counts},
    		warnings => []
	};
}

sub get_PATRIC_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["id","source","workspace"],{});
	my $data = Bio::KBase::ObjectAPI::utilities::rest_download({url => Bio::KBase::utilities::conf("GenomeImporter","data_api_url")."genome/?genome_id=".$args->{id}."&http_accept=application/json",token => Bio::KBase::utilities::token()});
	$data = $data->[0];
	$data = Bio::KBase::utilities::args($data,[],{
		genome_name => "Unknown",
		taxon_lineage_names => ["Unknown"],
		completion_date => undef
	});
	my $genome_input = {
		id => $args->{id},
		domain => $data->{taxon_lineage_names}->[0],
		source => "PATRIC",
		source_id => $args->{id},
		gc_content => $data->{gc_content},
		assembly_ref => undef,
		dna_size => 0,
		md5 => undef,
		scientific_name => $data->{genome_name},
		taxon_ref => undef,
		taxonomy => join(":",@{$data->{taxon_lineage_names}}),
		external_source_origination_date => $data->{completion_date},
		contig_ids => [],
		contig_lengths => [],
		num_contigs => 0,
		features => [],
		non_coding_features => [],
		feature_counts => {
	        CDS => 0,
	        gene => 0,
	        ncRNA => 0,
	        "non-protein_encoding_gene" => 0,
	        protein_encoding_gene => 0,
	        rRNA => 0,
	        regulatory => 0,
	        repeat_region => 0,
	        tRNA => 0,
	        tmRNA => 0
    		}
	};
	#Retrieving feature information
	my $start = 0;
	my $params = {};
	my $loopcount = 0;
	my $ftrcount = 0;
	my $allftrs = [];
	while ($start >= 0 && $loopcount < 100) {
		$loopcount++;#Insurance that no matter what, this loop won't run for more than 100 iterations
		my $ftrdata = Bio::KBase::ObjectAPI::utilities::rest_download({url => Bio::KBase::utilities::conf("GenomeImporter","data_api_url")."genome_feature/?genome_id=".$args->{id}."&http_accept=application/json&limit(10000,$start)",token => Bio::KBase::utilities::token()},$params);
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
	my $contig_feature_hash = {};
	my $aamd5 = [];
	my $namd5 = [];
	my $aamd5_feature_hash = {};
	my $namd5_feature_hash = {};
	for (my $i=0; $i < @{$allftrs}; $i++) {
		my $ftrdata = $allftrs->[$i];
		if ($ftrdata->{annotation} eq "PATRIC" 
			&& defined($ftrdata->{sequence_id}) 
			&& defined($ftrdata->{patric_id})
			&& defined($ftrdata->{feature_type})
			&& defined($ftrdata->{start})
			&& defined($ftrdata->{strand})
			&& defined($ftrdata->{na_sequence_md5})
			&& defined($ftrdata->{na_length})) {
			push(@{$namd5},$ftrdata->{na_sequence_md5});
			$namd5_feature_hash->{$ftrdata->{na_sequence_md5}} = $ftrdata;
			if (defined($ftrdata->{aa_sequence_md5})) {
				push(@{$aamd5},$ftrdata->{aa_sequence_md5});
				$aamd5_feature_hash->{$ftrdata->{aa_sequence_md5}} = $ftrdata;
			}
			push(@{$contig_feature_hash->{$ftrdata->{sequence_id}}},$ftrdata);	
		}
	}
	print "Feature contig count:".keys(%{$contig_feature_hash})."\n";
	#Pulling dna and protein sequences
	$self->query_for_sequences($namd5,$namd5_feature_hash,"na");
	$self->query_for_sequences($aamd5,$aamd5_feature_hash,"aa");
	#Pulling contig information
	if (!defined($d)) {
		$d = P3DataAPI->new();
	}
	my @res = $d->query("genome_sequence",
		["select", "genome_id", "accession", "sequence"],
		["eq", "genome_id", $args->{id}],
	);
	my $contigHash = {};
	for my $ent (@res) {
	    $contigHash->{$ent->{accession}} = $ent->{sequence};
	}
	my $contigarray;
	my $contigids = [];
	foreach my $contigid (keys(%{$contigHash})) {
		$genome_input->{num_contigs}++;
		push(@{$contigids},$contigid);
		push(@{$genome_input->{contig_lengths}},length($contigHash->{$contigid}));
		$genome_input->{dna_size} += length($contigHash->{$contigid});
		push(@{$contigarray},$contigHash->{$contigid});
	}
	$genome_input->{contig_ids} = [sort { $a cmp $b } @{$contigids}];
	print "Genome contig count:".$genome_input->{num_contigs}."\n";
	my $sortedcontigs = [sort { $a cmp $b } @{$contigarray}];
	my $str = "";
	for (my $i=0; $i < @{$sortedcontigs}; $i++) {
		if (length($str) > 0) {
			$str .= ";";
		}
		$str .= $sortedcontigs->[$i];	
	}
	$genome_input->{md5} = Digest::MD5::md5_hex($str);
	#Saving contigs
	$genome_input->{assembly_ref} = $args->{workspace}."/".$args->{id}.".contigs";
	$genome_input->{assembly_ref} = $self->save_assembly({
		workspace => $args->{workspace},
		data => $contigHash,
		name => $args->{id}.".contigs"
	});
	#Creating features array
	my $index = 1;
	my $added_contigs = {};
	foreach my $contig (@{$genome_input->{contig_ids}}) {
		if (defined($contig_feature_hash->{$contig})) {
			$added_contigs->{$contig} = 1;
			my $sortedftrlist = [sort { $a->{start} <=> $b->{start} } @{$contig_feature_hash->{$contig}}];
			for (my $i=0; $i < @{$sortedftrlist}; $i++) {
				if (defined($sortedftrlist->[$i]->{na_sequence})) {
					my $data = $sortedftrlist->[$i];
					my $ftrobj = {
						id => $genome_input->{id}."_".$index,
						type => $data->{feature_type},
						aliases=>[[$genome_input->{source}."_id",$data->{patric_id}]],
						location => [[$data->{sequence_id},$data->{start},$data->{strand},$data->{na_length}]],
					};
					$ftrobj->{location}->[0]->[1] += 0;
					if (defined($data->{product})) {
						$ftrobj->{functions} = [split(/\s*;\s+|\s+[\@\/]\s+/,$data->{product})];
						for (my $k=0; $k < @{$ftrobj->{functions}}; $k++) {
							$genome_input->{ontologies}->{SSO}->{$ftrobj->{id}}->{$ftrobj->{functions}->[$k]} = 1;
						}
					} else {
						$ftrobj->{functions} = ["Unknown"];
					}
					if (defined($data->{na_sequence})) {
						$ftrobj->{dna_sequence} = $data->{na_sequence};
						$ftrobj->{dna_sequence_length} = $data->{na_length};
						$ftrobj->{md5} = Digest::MD5::md5_hex($data->{na_sequence});
					}
					if (defined($data->{aa_sequence})) {
						$ftrobj->{type} = "gene";
						$ftrobj->{protein_translation} = $data->{aa_sequence};
						$ftrobj->{protein_translation_length} = $data->{aa_length};
						$ftrobj->{protein_md5} = $data->{aa_sequence_md5};
					}
					if (defined($data->{refseqgene}->{refseq_locus_tag})) {
						push(@{$ftrobj->{aliases}},["RefSeq",$data->{refseqgene}->{refseq_locus_tag}]);
					}
					if (defined($data->{refseqgene}->{gene})) {
						push(@{$ftrobj->{aliases}},["RefSeq",$data->{refseqgene}->{gene}]);
					}
					if (defined($data->{refseqgene}->{product})) {
						$genome_input->{ontologies}->{RefSeq}->{$ftrobj->{id}}->{$data->{refseqgene}->{product}} = 1;
					}
					if (defined($data->{figfam_id})) {
						$genome_input->{ontologies}->{FIGFAM}->{$ftrobj->{id}}->{$data->{figfam_id}} = 1;
					}
					if (defined($data->{pgfam_id})) {
						$genome_input->{ontologies}->{PGFAM}->{$ftrobj->{id}}->{$data->{pgfam_id}} = 1;
					}
					if (defined($data->{plfam_id})) {
						$genome_input->{ontologies}->{PLFAM}->{$ftrobj->{id}}->{$data->{plfam_id}} = 1;
					}
					if (defined($data->{go})) {
						for (my $k=0; $k < @{$data->{go}}; $k++) {
							my $array = [split(/\|/,$data->{go}->[$k])];
							$genome_input->{ontologies}->{GO}->{$ftrobj->{id}}->{$array->[0]} = 1;
						}
					}
					$index++;
					if ($ftrobj->{type} eq "gene") {
						push(@{$genome_input->{features}},$ftrobj);
						$genome_input->{feature_counts}->{CDS} += 1;
						$genome_input->{feature_counts}->{gene} += 1;
					} else {
						push(@{$genome_input->{non_coding_features}},$ftrobj);
						if ($ftrobj->{type} eq "tRNA") {
							$genome_input->{feature_counts}->{tRNA} += 1;
						} elsif ($ftrobj->{type} eq "rRNA") {
							$genome_input->{feature_counts}->{rRNA} += 1;
						} elsif ($ftrobj->{type} eq "rRNA") {
							$genome_input->{feature_counts}->{"non-protein_encoding_gene"} += 1;
						}
					}
				} else {
					print $sortedftrlist->[$i]->{patric_id}." no DNA sequence!\n";
				}
			}
		} else {
			print $contig." has no features!\n";
		}
	}
	foreach my $contig (keys(%{$contig_feature_hash})) {
		if (!defined($added_contigs->{$contig})) {
			my $ftrcount = @{$contig_feature_hash->{$contig}};
			print $contig." ".$ftrcount." features never added!\n";
		}
	}
	#Saving genome
	$self->save_genome({
		workspace => $args->{workspace},
		data => $genome_input,
		name => $genome_input->{id}
	});
}

sub get_SEED_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["id","source","workspace"],{});
	#Configuring the server based on the requested source
	my $sapsvr;
	if ($args->{source} eq "pubseed") {
		$sapsvr = Bio::ModelSEED::Client::SAP->new();
	} elsif ($args->{source} eq "coreseed") {
		$sapsvr = Bio::ModelSEED::Client::SAP->new({url => "https://core.theseed.org/FIG/sap_server.cgi"});
	}
	my $translation = {
		pubseed => "PubSEED",
		coreseed => "CoreSEED"
	};
	my $data = $sapsvr->genome_data({
		-ids => [$args->{id}],
		-data => [qw(gc-content dna-size name taxonomy domain genetic-code)]
	});
	if (!defined($data->{$args->{id}})) {
    		$self->error("PubSEED genome ".$args->{id}." not found!");
    }
    #Initializing the genome input
    my $genome_input = {
		id => $args->{id},
		domain => $data->{$args->{id}}->[4],
		source => $translation->{$args->{source}},
		source_id => $args->{id},
		gc_content => $data->{$args->{id}}->[0]/100,
		assembly_ref => undef,
		dna_size => $data->{$args->{id}}->[1]+0,
		md5 => undef,
		scientific_name => $data->{$args->{id}}->[2],
		genetic_code => $data->{$args->{id}}->[5]+0,
		taxon_ref => undef,
		taxonomy => $data->{$args->{id}}->[3],
		external_source_origination_date => $data->{completion_date},
		num_contigs => 0,
		contig_lengths => [],
		contig_ids => [],
		features => [],
		non_coding_features => [],
		feature_counts => {
	        CDS => 0,
	        gene => 0,
	        ncRNA => 0,
	        "non-protein_encoding_gene" => 0,
	        protein_encoding_gene => 0,
	        rRNA => 0,
	        regulatory => 0,
	        repeat_region => 0,
	        tRNA => 0,
	        tmRNA => 0
    		}
	};
    #Processing the assembly
    my $genomeHash = $sapsvr->genome_contigs({
		-ids => [$args->{id}]
	});
    my $contigList = $genomeHash->{$args->{id}};
    my $contigHash = $sapsvr->contig_sequences({
		-ids => $contigList
	});
    my $contigarray;
    my $newhash = {};
    my $contigids = [];
	foreach my $contigid (keys(%{$contigHash})) {
		my $array = [split(/:/,$contigid)];
		$newhash->{$array->[1]} = $contigHash->{$contigid};
		push(@{$contigarray},$contigHash->{$contigid});
		$genome_input->{num_contigs}++;
		push(@{$genome_input->{contig_lengths}},length($contigHash->{$contigid}));
		push(@{$contigids},$array->[1]);
	}
	$genome_input->{contig_ids} = [sort { $a cmp $b } @{$contigids}];
	my $sortedcontigs = [sort { $a cmp $b } @{$contigarray}];
	my $str = "";
	for (my $i=0; $i < @{$sortedcontigs}; $i++) {
		if (length($str) > 0) {
			$str .= ";";
		}
		$str .= $sortedcontigs->[$i];	
	}
	$genome_input->{md5} = Digest::MD5::md5_hex($str);
    	#Saving contigs
    	$genome_input->{assembly_ref} = $args->{workspace}."/".$args->{id}.".contigs";
	$genome_input->{assembly_ref} = $self->save_assembly({
		workspace => $args->{workspace},
		data => $newhash,
		name => $args->{id}.".contigs"
	});
	#Pulling feature data
	my $featureHash = $sapsvr->all_features({-ids => $args->{id}});
	my $featureList = $featureHash->{$args->{id}};
	my $functions = $sapsvr->ids_to_functions({-ids => $featureList});
	my $dnaseq = $sapsvr->ids_to_sequences({-ids => $featureList});
	my $locations = $sapsvr->fid_locations({-ids => $featureList});
	my $sequences = $sapsvr->fids_to_proteins({-ids => $featureList,-sequence => 1});
	my $contig_feature_hash = {};
	for (my $i=0; $i < @{$featureList}; $i++) {
		my $feature = {
			type => "unknown",
			aliases => [[$translation->{$args->{source}}."_id",$featureList->[$i]]],
  		};
  		if ($featureList->[$i] =~ m/\.([^\.]+)\.\d+$/) {
  			$feature->{type} = $1;
  		}
  		if (defined($functions->{$featureList->[$i]})) {
  			$feature->{functions} = [split(/\s*;\s+|\s+[\@\/]\s+/,$functions->{$featureList->[$i]})];
		}
		if (defined($sequences->{$featureList->[$i]})) {
			$feature->{type} = "gene";
			$feature->{protein_translation} = $sequences->{$featureList->[$i]};
			$feature->{protein_translation_length} = length($feature->{protein_translation});
  			$feature->{dna_sequence_length} = 3*$feature->{protein_translation_length};
  			$feature->{protein_md5} = Digest::MD5::md5_hex($feature->{protein_translation});
		}
		if (defined($dnaseq->{$featureList->[$i]})) {
			$feature->{dna_sequence} = $dnaseq->{$featureList->[$i]};
			$feature->{md5} = Digest::MD5::md5_hex($feature->{dna_sequence});
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
		#print $feature->{type}."\n";
		#print $feature->{dna_sequence}."\n";
		if (defined($feature->{location})
			&& defined($feature->{dna_sequence})
			&& defined($feature->{type})
			&& $feature->{type} ne "opr"
			&& $feature->{type} ne "pbs"
			&& $feature->{type} ne "prm"
			&& $feature->{type} ne "trm") {
  			push(@{$contig_feature_hash->{$feature->{location}->[0]->[0]}},$feature);
		}
	}
	#Creating features array
	my $index = 1;
	foreach my $contig (@{$genome_input->{contig_ids}}) {
		print "Contig:".$contig."\n";
		if (defined($contig_feature_hash->{$contig})) {
			my $sortedftrlist = [sort { $a->{location}->[0]->[1] <=> $b->{location}->[0]->[1] } @{$contig_feature_hash->{$contig}} ];
			for (my $i=0; $i < @{$sortedftrlist}; $i++) {
				my $ftrobj = $sortedftrlist->[$i];
				$ftrobj->{id} = $genome_input->{id}."_".$index;
				if (defined($ftrobj->{functions})) {
		  			for (my $k=0; $k < @{$ftrobj->{functions}}; $k++) {
						$genome_input->{ontologies}->{SSO}->{$ftrobj->{id}}->{$ftrobj->{functions}->[$k]} = 1;
					}
				} else {
					$ftrobj->{functions} = ["Unknown"];
				}
				$index++;
				if ($ftrobj->{type} eq "gene") {
					push(@{$genome_input->{features}},$ftrobj);
					$genome_input->{feature_counts}->{CDS} += 1;
					$genome_input->{feature_counts}->{gene} += 1;
				} else {
					push(@{$genome_input->{non_coding_features}},$ftrobj);
					if ($ftrobj->{type} eq "rna") {
						$genome_input->{feature_counts}->{rRNA} += 1;
					} else {
						$genome_input->{feature_counts}->{"non-protein_encoding_gene"} += 1;
					}
				}
			}
		}
	}
	#Saving genome
	$self->save_genome({
		workspace => $args->{workspace},
		data => $genome_input,
		name => $genome_input->{id}
	});
}

sub save_genome {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["workspace","data"],{
		name => $args->{data}->{id}
	});
	#Checking if the id is a relevant tax id
	my $taxid = $args->{data}->{id};
	$taxid =~ s/\.\d+//;
	eval {
		Bio::KBase::kbaseenv::get_object_info([{"ref" => "ReferenceTaxons/".$taxid."_taxon"}],0);
	};
	if (!$@) {
		$args->{data}->{taxon_ref} = 	"ReferenceTaxons/".$taxid."_taxon";
	}
	#Building genome
	my $genome = $self->get_base_genome($args->{data});
	#Automatically creating CDS for all coding genes
	foreach my $ftr (@{$genome->{features}}) {
		my $cdsftr = {};
		foreach my $field (keys(%{$ftr})) {
			$cdsftr->{$field} = $ftr->{$field};
		}
		$cdsftr->{id} = $ftr->{id}."_CDS_1";
		$cdsftr->{type} = "CDS";
		$cdsftr->{parent_gene} = $ftr->{id};
		$ftr->{cdss} = [$cdsftr->{id}];
		push(@{$genome->{cdss}},$cdsftr);
	}
	#Building ontology events
	my $input = {
		object => $genome,
		events => [],
		save => 1,
		output_name => $args->{data}->{id},
		output_workspace => $args->{workspace},
		type => "KBaseGenomes.Genome"
	};
	foreach my $ontology (keys(%{$args->{data}->{ontologies}})) {
		my $newevent = {
			description => $ontology." annotations imported from ".$args->{data}->{source},
			ontology_id => $ontology,
			method => "GenomeImporter-import_external_genome",
			method_version => "1.0",
			timestamp => Bio::KBase::utilities::timestamp(1),
			ontology_terms => {}
		};
		foreach my $gene (keys(%{$args->{data}->{ontologies}->{$ontology}})) {
			foreach my $term (keys(%{$args->{data}->{ontologies}->{$ontology}->{$gene}})) {
				push(@{$newevent->{ontology_terms}->{$gene}},{term => $term});
			} 
		}
		push(@{$input->{events}},$newevent);
	}
	Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/ontology_api_input.json",[Bio::KBase::utilities::to_json($input,1)]);
	#Saving the genome using the ontology service
	print "Calling ontology service!";
	my $anno_ontology_client = annotation_ontology_api::annotation_ontology_apiServiceClient->new(undef,token => Bio::KBase::utilities::token());
	my $output = $anno_ontology_client->add_annotation_ontology_events($input);
	Bio::KBase::kbaseenv::add_object_created({
		"ref" => $output->{output_ref},
		"description" => "Genome object for ".$args->{name}
	});
	return $output->{output_ref};
}

sub save_assembly {
	my($self,$args) = @_;
	$args = Bio::KBase::utilities::args($args,["workspace","data","name"],{});
	my $filename = Bio::KBase::utilities::conf("GenomeImporter","scratch")."/".$args->{name}.".fa";
	open(my $fo, ">", $filename) || die "Could not open file ".$filename;
	foreach my $contigid (keys(%{$args->{data}})) {
		print $fo ">".$contigid."\n";
		my $sequence = $args->{data}->{$contigid};
		my $seq_len = length($sequence);
		for (my $i = 0; $i < $seq_len; $i += 60) {
		    my $segment = substr($sequence, $i, 60);
		    print $fo $segment."\n";
		}
	}
	close($fo);
	my $ga = Bio::KBase::kbaseenv::ac_client();
	my $assemblyref = $ga->save_assembly_from_fasta({
		file => {path => $filename},
        workspace_name => $args->{workspace},
        assembly_name => $args->{name}
	});
	Bio::KBase::kbaseenv::add_object_created({
		"ref" => $assemblyref,
		"description" => "Contigs for ".$args->{name}
	});
	return $assemblyref;
}

sub util_initialize_call {
	my ($self,$params,$ctx) = @_;
	print "Import parameters:".Bio::KBase::ObjectAPI::utilities::TOJSON($params,1);
	if (defined($ctx)) {
		Bio::KBase::kbaseenv::initialize_call($ctx);
	}
	return $params;
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
	$self->util_initialize_call($params,$ctx);
	$args = Bio::KBase::utilities::args($args,["genome_ids","workspace","source"],{});
    my $genomes = [split(/[\n;\|]+/,$args->{genome_ids})];
    my $htmlmessage = "<p>";
    for (my $i=0; $i<@{$genomes};$i++) {
    	print "Now importing ".$genomes->[$i]." from ".$args->{source}."\n";
    	eval {
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
    };
    	if ($@) {
			$htmlmessage .= $genomes->[$i]." failed!<br>".$@;
		} else {
			$htmlmessage .= $genomes->[$i]." succeeded!<br>";
		}
    }
    $htmlmessage .= "</p>";
	Bio::KBase::utilities::print_report_message({
		message => $htmlmessage,html=>1,append => 0
	});
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
