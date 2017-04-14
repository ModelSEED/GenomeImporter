use GenomeImporter::GenomeImporterImpl;

use GenomeImporter::GenomeImporterServer;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = GenomeImporter::GenomeImporterImpl->new;
    push(@dispatch, 'GenomeImporter' => $obj);
}


my $server = GenomeImporter::GenomeImporterServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
