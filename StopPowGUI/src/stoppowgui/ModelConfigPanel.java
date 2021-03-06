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

import javax.swing.JDialog;

/**
 * Generic class for model configuration panels
 * @author Alex Zylstra
 * @date 2013/06/05
 */
public abstract class ModelConfigPanel extends javax.swing.JPanel {
    /** The ModelManager in use by this application */
    ModelManager models;
    /** The parent dialog for this panel */
    JDialog dialog;
    /** Whether the model is already created. If this is set to true,
        the expected value for extending classes is that they will
        update an existing model. If it is false (default) then a new
        model is created. */
    boolean editMode;
    /** Previous key used for this model. Used with the edit functionality
     * to prevent accidentally overwriting another unrelated model
     */
    String previousKey;
    
    /**
     * Construct a new generic model configuration panel.
     * @param models the ModelManager to be used
     */
    public ModelConfigPanel(ModelManager models) {
        this.models = models;   
        editMode = false;
        previousKey = "";
    }
    
    /**
     * Set a parent dialog for this panel.
     * This is necessary for certain buttons to work properly
     * @param dialog the Dialog that this panel is displayed in
     */
    public void setDialog(JDialog dialog)
    {
        this.dialog = dialog;
    }
    
    /** 
     * Get the type of model created by this panel
     * @return The type of model created by this panel
     */
    abstract public String get_type();
    
    /**
     * Get if the panel is in edit mode.
     * @return True if it is in edit mode.
     */
    public boolean isEditMode()
    {
        return editMode;
    }
    
    /**
     * Get the previous key used for the model configured
     * by this panel.
     * @return The key string.
     */
    public String getPreviousKey()
    {
        return previousKey;
    }
}
