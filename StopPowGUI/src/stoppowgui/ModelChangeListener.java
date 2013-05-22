package stoppowgui;

/**
 * Implement a listener for model changes.
 * If a stopping power model in the ModelManger is changed,
 * added, or removed, the model_changed function is fired
 * and given the key (i.e. name) of the affected model.
 * 
 * @brief A Listener for dE/dx model changes
 * @class ModelChangeListener
 * @date 2013/05/10
 * @author Alex Zylstra
 */
public interface ModelChangeListener extends java.util.EventListener
{
    /**
     * Function called when a model is changed
     * @param evt A ModelChangeEvent corresponding to this event
     */
    void model_changed(ModelChangeEvent evt);
    
    /**
     * Function called when a model is added
     * @param evt A ModelChangeEvent corresponding to this event
     */
    void model_added(ModelChangeEvent evt);
    
    /**
     * Function called when a model is deleted
     * @param evt A ModelChangeEvent corresponding to this event
     */
    void model_deleted(ModelChangeEvent evt);
}
