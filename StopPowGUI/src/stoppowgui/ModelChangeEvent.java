package stoppowgui;

/**
 * An Event for notifying Listeners that a stopping
 * power model has been changed.
 * @author alex
 */
public class ModelChangeEvent extends java.util.EventObject
{
    /** key String for the model that was changed */
    String key;
    
    /**
     * Create a new event.
     * @param source The originating source of the event
     * @param key_in The Key String for the originating model
     */
    ModelChangeEvent(Object source, String key_in)
    {
        super(source);
        key = key_in;
    }
    
    /**
     * Get the key for the originating event
     * @return key String
     */
    String get_key()
    {
        return key;
    }
}

