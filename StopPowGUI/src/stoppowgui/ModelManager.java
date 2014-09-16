// StopPow - a charged-particle stopping power library
// Copyright (C) 2014  Massachusetts Institute of Technology / Alex Zylstra

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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
 * @date 2013/06/05
 */
public class ModelManager 
{
    /** Use a Map data structure for storing the data */
    private Map<String,StopPow> models;
    /** Also store configuration dialog panels for each model */
    private Map<String,ModelConfigPanel> modelConfig;
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
        models = new HashMap<String,StopPow>();
        modelConfig = new HashMap<String,ModelConfigPanel>();
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
            for( Map.Entry<String,StopPow> row : models.entrySet() )
            {
                // have to get the "value" from the Map, then the model
                // from it (since each value is a ModelEntry object)
                ret.add( row.getValue() );
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
     * Add a model to the ModelManager
     * @param name The name of the model (will be used as key - must be unique!)
     * @param new_model the model (a StopPow object)
     * @param configPanel the configuration panel dialog for changing the model
     */
    public void add_model(String name, StopPow new_model, ModelConfigPanel configPanel) 
            throws java.lang.IllegalArgumentException
    {
        // sanity check:
        if( models.containsKey(name) )
            throw new java.lang.IllegalArgumentException("Key already exists");
        
        // acquire the lock for the manager only:
        lockManager();
        try
        {
            models.put(name , new_model);
            modelConfig.put(name, configPanel);
            modelLocks.put(name, new ReentrantLock());
        }
        finally
        {
            unlockManager();
        }
        
        // fire change listeners:
        this.fireModelAdded(name);
    }
    
    public void change_model(String name, StopPow new_model, ModelConfigPanel new_panel)
            throws java.lang.IllegalArgumentException
    {
        // sanity check:
        if( !models.containsKey(name) )
            throw new java.lang.IllegalArgumentException("Key does not exist!");
        
        // acquire the lock for the manager and this key:
        lockManager();
        lockModel(name);
        try
        {
            models.put(name, new_model);
            modelConfig.put(name, new_panel);
        }
        finally
        {
            unlockManager();
            unlockModel(name);
        }
        
        // fire change listeners:
        this.fireModelChange(name);
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
                modelConfig.remove(key);
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
        return models.get(key);
    }
    
    /**
     * Get the configuration panel for the model corresponding to the given key
     * @param key The key value
     * @return The configuration panel object
     */
    public ModelConfigPanel get_panel(String key)
    {
        return modelConfig.get(key);
    }
    
    /**
     * Get the type of model corresponding to a given key
     * @param key The key (i.e. name) of model you want
     * @return The type of model
     */
    public String get_type(String key)
    {
        return models.get(key).get_type();
    }
    
    /**
     * Get the info for model corresponding to a given key
     * @param key The key (i.e. name) of model you want
     * @return The info string for that model
     */
    public String get_info(String key)
    {
        return models.get(key).get_info();
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

