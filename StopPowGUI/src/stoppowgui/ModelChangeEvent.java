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

