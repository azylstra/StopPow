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

package stoppowgui.util;

import java.awt.Window;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import javax.swing.JToggleButton;

/**
 *
 * @author alex
 */
public class ToggleButtonWindowListener implements WindowListener {
    /** The toggle button associated with display of the window */
    JToggleButton button;
    
    /**
     * Create a new ToggleButtonWindowListener.
     * @param button The button to associate with the state of the window
     */
    public ToggleButtonWindowListener(JToggleButton button) {
        this.button = button;
    }
    
    @Override
    public void windowOpened(WindowEvent e) {
        // if the window is opened via other means,
        // set button state to true:
        button.setSelected(true);
    }

    @Override
    public void windowClosing(WindowEvent e) {
        // if the window is closed via other means, make sure button state
        // is set to false:
        button.setSelected(false);
    }

    @Override
    public void windowClosed(WindowEvent e) {
        // if the window is closed via other means, make sure button state
        // is set to false:
        button.setSelected(false);
    }

    @Override
    public void windowIconified(WindowEvent e) {
        // do nothing
    }

    @Override
    public void windowDeiconified(WindowEvent e) {
        // do nothing
    }

    @Override
    public void windowActivated(WindowEvent e) {
        // do nothing
    }

    @Override
    public void windowDeactivated(WindowEvent e) {
        // do nothing
    }
    
}
