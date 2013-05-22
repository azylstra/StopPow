package stoppowgui;

import cStopPow.StopPow;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 * This class will control all stopping power models used in the code.
 * It implements functionality for thread-safe additions and deletions
 * of models, plus functionality for firing Listener events, so that
 * other code can be notified when the models are changed.
 * 
 * @class ModelManager
 * @author Alex Zylstra
 * @date 2013/05/18
 */
public class ModelManager 
{
    /** Use a Map data structure for storing the data */
    private Map<String,ModelEntry> models;
    /** Provide mutex locking functionality for individual models */
    private Map<String,Lock> modelLocks;
    /** For mutex locking the whole ModelManager */
    private final Lock lock = new ReentrantLock();
    /** A list of Listeners for change notifications */
    ArrayList<ModelChangeListener> listeners;

    /** Constructor */
    public ModelManager()
    {
        // initialize:
        models = new HashMap<String,ModelEntry>();
        modelLocks = new HashMap<String,Lock>();
        listeners = new ArrayList<ModelChangeListener>();
    }
    
    // ---------------------------------------
    //      Access or modify model data
    // ---------------------------------------
    /**
     * Get all of the models.
     * @return an ArrayList representation of all models
     */
    public ArrayList<StopPow> get_models()
    {
        // build an ArrayList to return:
        ArrayList<StopPow> ret;
        
        // acquire the lock for the manager
        lockManager();
        // do work in try/finally to ensure that lock is released
        // at the end:
        try
        {
            ret = new ArrayList<StopPow>();

            // iterate over all values in the Map and add to the ArrayList:
            for( Map.Entry<String,ModelEntry> row : models.entrySet() )
            {
                // have to get the "value" from the Map, then the model
                // from it (since each value is a ModelEntry object)
                ret.add( row.getValue().model );
            }
        }
        finally
        {        
            // release the lock:
            unlockManager();
        }
        
        return ret;
    }        
    
    /**
     * Add a model to the ModelManager.
     * @param name The name of the model (will be used as key - must be unique!)
     * @param model_name The type of model, i.e. SRIM, Li-Petrasso, ...
     * @param model The stopping power model itself
     * @param info Any additional information about the model
     */
    public void add_model(String name, String model_name, StopPow model, String info)
            throws java.lang.IllegalArgumentException
    {        
        ModelEntry new_entry = new ModelEntry(name,model_name,model,info);
        add_model(new_entry);
    }
    
    /**
     * Add a model to the ModelManager
     * @param new_model the model, wrapped as ModelEntry
     */
    public void add_model(ModelEntry new_model) throws java.lang.IllegalArgumentException
    {
        // sanity check:
        if( models.containsKey(new_model.name) )
            throw new java.lang.IllegalArgumentException("Key already exists");
        
        // acquire the lock for the manager only:
        lockManager();
        try
        {
            models.put(new_model.name , new_model);
        }
        finally
        {
            unlockManager();
        }
        
        // fire change listeners:
        this.fireModelAdded(new_model.name);
    }
    
    /**
     * Remove a model from the ModelManager
     * @param key the key value corresponding to the model you wish to remove
     */
    public void remove_model(String key)
    {
        // sanity check:
        if( !models.containsKey(key) )
            return;
        
        // acquire the lock for the manager:
        lockManager();
        try
        {
            // also lock the model:
            lockModel(key);
            try {
                models.remove(key);
            }
            finally {
                unlockModel(key);
                // after unlocking, remove the lock:
                modelLocks.remove(key);
            }
        }
        finally
        {
            unlockManager();
        }
        
        // fire change listeners:
        this.fireModelDeleted(key);
    }
        
    /**
     * Get the model corresponding to a given key
     * @param key The key (i.e. name) of model you want
     * @return The StopPow object corresponding to key
     */
    public StopPow get_model(String key)
    {
        return models.get(key).model;
    }
    
    /**
     * Get the model as a ModelEntry
     * @param key the Key (i.e. name) of the model you want
     * @return the ModelEntry object corresponding to that model
     */
    public ModelEntry get_model_entry(String key)
    {
        return models.get(key);
    }
    
    /**
     * Get the type of model corresponding to a given key
     * @param key The key (i.e. name) of model you want
     * @return The type of model
     */
    public String get_type(String key)
    {
        return models.get(key).type;
    }
    
    /**
     * Get the info for model corresponding to a given key
     * @param key The key (i.e. name) of model you want
     * @return The info string for that model
     */
    public String get_info(String key)
    {
        return models.get(key).info;
    }
        
    /**
     * Check if a key containsKey in the ModelManager already
     * @param key The String you wish to check
     * @return true if key is already used, false otherwise
     */
    public boolean containsKey(String key)
    {
        return models.containsKey(key);
    }
    
    /**
     * Get the Set of all keys being used.
     * @return Set<String> of all keys
     */
    public Set<String> get_keys()
    {
        return models.keySet();
    }
    
    // ---------------------------------------
    //       Mutex Lock Functionality
    // ---------------------------------------
    /**
     * Acquire mutex lock for the model manager
     */
    public void lockManager()
    {
        lock.lock();
    }
    /**
     * Try to acquire mutex lock for the model manager
     * @return true if lock was acquired
     */
    public boolean tryLockManager()
    {
        return lock.tryLock();
    }
    /**
     * Release lock for the model manager
     */
    public void unlockManager()
    {
        lock.unlock();
    }
    
    /**
     * Acquire a mutex lock for a specific model. If the given key
     * does not correspond to any models in the ModelManager,
     * this method does nothing.
     * @param key the String key for the model to lock
     */
    public void lockModel(String key)
    {
        // sanity check:
        if( modelLocks.containsKey(key) )
        {
            Lock l = modelLocks.get(key);
            l.lock();
        }
    }
    /**
     * Release mutex lock for a specific model. If the given key
     * does not correspond to any models in the ModelManager,
     * this method does nothing.
     * @param key the String key for the model to unlock
     */
    public void unlockModel(String key)
    {
        // sanity check:
        if( modelLocks.containsKey(key) )
        {
            Lock l = modelLocks.get(key);
            l.unlock();
        }
    }
    /**
     * Acquire mutex lock for all models.
     */
    public void lockAllModels()
    {
        for( Lock iterLock : modelLocks.values() )
            iterLock.lock();
    }
    /**
     * Release mutex lock for all models.
     */
    public void unlockAllModels()
    {
        for( Lock iterLock : modelLocks.values() )
            iterLock.unlock();
    }
    
    // ---------------------------------------
    //       Listeners Functionality
    // ---------------------------------------
    /**
     * Add a listener to the ModelManager. It will get called whenver
     * any models are changed, added, or removed.
     * @param listener 
     */
    public void addModelChangeListener(ModelChangeListener listener)
    {
        listeners.add(listener);
    }
    
    /**
     * Fire model changed events to listeners
     * @param key Key of the model changed
     */
    private void fireModelChange(String key)
    {
        // construct an event:
        ModelChangeEvent evt = new ModelChangeEvent(this,key);
        // loop over all listeners and notify them:
        for( ModelChangeListener l : listeners )
            l.model_changed(evt);
    }
    
    /**
     * Fire model added events to listeners
     * @param key Key of the model added
     */
    private void fireModelAdded(String key)
    {
        // construct an event:
        ModelChangeEvent evt = new ModelChangeEvent(this,key);
        // loop over all listeners and notify them:
        for( ModelChangeListener l : listeners )
            l.model_added(evt);
    }
    
    /**
     * Fire model deleted events to listeners
     * @param key Key of the model deleted
     */
    private void fireModelDeleted(String key)
    {
        // construct an event:
        ModelChangeEvent evt = new ModelChangeEvent(this,key);
        // loop over all listeners and notify them:
        for( ModelChangeListener l : listeners )
            l.model_deleted(evt);
    }
}

