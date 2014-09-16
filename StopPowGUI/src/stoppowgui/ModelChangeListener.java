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
