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
