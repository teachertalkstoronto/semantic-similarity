//import java.util.LinkedHashSet;
import java.util.Set;
import org.openrdf.model.URI;
import slib.sml.sm.core.utils.SMconf;

public class Gene {
    public int id;
    public Set<URI> annot;
    public Gene(int i, Set<URI> a) {
        id  = i;
        annot = a;
    }

}