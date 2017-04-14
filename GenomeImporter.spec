/*
A KBase module: GenomeImporter
This sample module contains one small method - filter_contigs.
*/

module GenomeImporter {
    /*
        Input parameters for the import_external_genome function
    */
    typedef structure {
        string genome_ids;
        string source;
        string workspace;
    } ImportGenomeParams;


    /*
        Output structure for the import_external_genome function
    */
    typedef structure {
        string report_name;
        string report_ref;
    } ImportGenomeResults;
    
    /*
        Function to import a list of genomes from a specified source
    */
    funcdef import_external_genome(ImportGenomeParams params)
        returns (ImportGenomeResults output) authentication required;
};
