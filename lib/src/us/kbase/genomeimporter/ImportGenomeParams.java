
package us.kbase.genomeimporter;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: ImportGenomeParams</p>
 * <pre>
 * Input parameters for the import_external_genome function
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "genome_ids",
    "source",
    "workspace"
})
public class ImportGenomeParams {

    @JsonProperty("genome_ids")
    private String genomeIds;
    @JsonProperty("source")
    private String source;
    @JsonProperty("workspace")
    private String workspace;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("genome_ids")
    public String getGenomeIds() {
        return genomeIds;
    }

    @JsonProperty("genome_ids")
    public void setGenomeIds(String genomeIds) {
        this.genomeIds = genomeIds;
    }

    public ImportGenomeParams withGenomeIds(String genomeIds) {
        this.genomeIds = genomeIds;
        return this;
    }

    @JsonProperty("source")
    public String getSource() {
        return source;
    }

    @JsonProperty("source")
    public void setSource(String source) {
        this.source = source;
    }

    public ImportGenomeParams withSource(String source) {
        this.source = source;
        return this;
    }

    @JsonProperty("workspace")
    public String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(String workspace) {
        this.workspace = workspace;
    }

    public ImportGenomeParams withWorkspace(String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((("ImportGenomeParams"+" [genomeIds=")+ genomeIds)+", source=")+ source)+", workspace=")+ workspace)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
