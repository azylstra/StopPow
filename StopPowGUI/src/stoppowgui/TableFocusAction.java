package stoppowgui;

import java.awt.Component;
import java.awt.event.*;
import java.beans.PropertyChangeListener;
import javax.swing.*;

/**
 * Provide a wrapper for javax.swing.Action, which is used to keep
 * track of the originating action and thisComponent.
 * @author alex
 */
abstract class WrappedAction implements Action
{
	private Action defaultAction;
	private JComponent thisComponent;
	private Object actionBinding;

        /**
         * Construct a new WrappedAction.
         * @param thisComponent the originating thisComponent
         * @param keyStroke the originating keystroke
         */
	public WrappedAction(JComponent component, KeyStroke keyStroke)
	{
            // set the class variable thisComponent:
            this.thisComponent = component;
            
            // local temp copy:
            Object tempKey = getKeyForActionMap(component, keyStroke);
            
            // sanity check:
            if (tempKey == null)
            {
                    String message = "No input mapping in WrappedAction for KeyStroke: " + keyStroke;
                    throw new IllegalArgumentException(message);
            }
            
            // actually set the key:
            setActionForKey( tempKey );
	}

        /**
         * Construct a new WrappedAction from thisComponent and key
         * @param thisComponent the originating thisComponent
         * @param actionBinding the key
         */
	public WrappedAction(JComponent component, Object actionKey)
	{
            this.thisComponent = component;
            setActionForKey( actionKey );
	}

        /**
         * Convert from KeyStroke to binding
         * @param thisComponent The originating thisComponent
         * @param keyStroke The key stroke triggering this event
         * @return 
         */
	private Object getKeyForActionMap(JComponent component, KeyStroke keyStroke)
	{
            // loop over three input maps:
            for (int i = 0; i < 3; i++)
            {
                InputMap inputMap = component.getInputMap(i);

                    if (inputMap != null) // sanity
                    {
                        // get the binding:
                        Object key = inputMap.get(keyStroke);

                        if (key != null)
                                return key;
                    }
            }

		return null;
	}

        /**
         * Replace a current binding for a key with a custom action
         * as defined in this class.
         * @param actionBinding the current binding
         */
	private void setActionForKey(Object actionKey)
	{
            //  Save the original Action

            this.actionBinding = actionKey;
            defaultAction = thisComponent.getActionMap().get(actionKey);

            if (defaultAction == null)
            {
                    String message = "no Action for action key: " + actionKey;
                    throw new IllegalArgumentException(message);
            }

            //  Replace the existing Action with this class
            install();
	}

        /**
         * Allow another class to invoke the original or default action.
         * @param e the triggering ActionEvent we want to use to invoke the default Action
         */
	public void invokeOriginalAction(ActionEvent e)
	{
            defaultAction.actionPerformed(e);
	}

        /**
         * Replace the default action with the custom action defined.
         */
	public void install()
	{
            thisComponent.getActionMap().put(actionBinding, this);
	}

	/**
         * Remove the custom action and re-enable the default.
         */
	public void unInstall()
	{
            thisComponent.getActionMap().put(actionBinding, defaultAction);
	}
//
//  Delegate the Action interface methods to the original Action
//
        /**
         * Add a PropertyChangeListener. This is applied to the default Action.
         * @param listener the listener you want to add
         */
	public void addPropertyChangeListener(PropertyChangeListener listener)
	{
		defaultAction.addPropertyChangeListener(listener);
	}

        /**
         * Get a property from this object using the key. This method
         * invokes the default action.
         * @param key The key value to query
         * @return The resulting property
         */
        @Override
	public Object getValue(String key)
	{
		return defaultAction.getValue(key);
	}

        /**
         * Check if the default action is enabled.
         * @return true if the default action is enabled, false otherwise
         */
        @Override
	public boolean isEnabled()
	{
		return defaultAction.isEnabled();
	}

        /**
         * Add a new property. This is applied to the default Action.
         * @param key a new key
         * @param newValue a new value
         */
        @Override
	public void putValue(String key, Object newValue)
	{
		defaultAction.putValue(key, newValue);
	}

        /**
         * Remove a specific PropertyChangeListener. Applied to the default action.
         * @param listener The listener to remove
         */
        @Override
	public void removePropertyChangeListener(PropertyChangeListener listener)
	{
		defaultAction.removePropertyChangeListener(listener);
	}

        /**
         * Toggle whether the default action should be used.
         * @param newValue set true to enable the default action.
         */
        @Override
	public void setEnabled(boolean newValue)
	{
		defaultAction.setEnabled(newValue);
	}

        /**
         * Get all keys associated with the default action.
         * @return all keys as an Object[]
         */
	public Object[] getKeys()
	{
		if (defaultAction instanceof AbstractAction)
		{
			AbstractAction abstractAction = (AbstractAction)defaultAction;
			return abstractAction.getKeys();
		}

		return null;
	}

        /**
         * Get a list of all property listeners
         * @return All listeners attached to the default action, null if
         * the default action is not an instance of AbstractActions
         */
	public PropertyChangeListener[] getPropertyChangeListeners()
	{
		if (defaultAction instanceof AbstractAction)
		{
			AbstractAction abstractAction = (AbstractAction)defaultAction;
			return abstractAction.getPropertyChangeListeners();
		}

		return null;
	}
}
/**
 * Use the above WrappedAction functionality to extend key bindings
 * for a JTable. Primarily, this supports letting a user tab through
 * a JTable, and then tab navigate to another component when the lower
 * right hand side of a table is reached.
 * @author Alex Zylstra
 * @date 2013/05/16
 */
public class TableFocusAction extends WrappedAction implements ActionListener
{
	private JTable table;
        private Component next_comp;

        /**
         * Construct a new TableFocusAction. The user must specify the table, the
         * keystroke, and which component to move to next.
         * @param table
         * @param keyStroke
         * @param next 
         */
	public TableFocusAction(JTable table, KeyStroke keyStroke, Component next)
	{
		super(table, keyStroke);
		this.table = table;
                this.next_comp = next;
	}

	/**
         * Implement the custom action.
         * @param e The ActionEvent that triggered this call.
         */
	public void actionPerformed(ActionEvent e)
	{
            // get some info on the dimensions of the table and where we are:
            int originalRow = table.getSelectedRow();
            int originalColumn = table.getSelectedColumn();
            int num_rows = table.getRowCount();
            int num_cols = table.getColumnCount();
            
            // current location in the table:
            int row = table.getSelectedRow();
            int col = table.getSelectedColumn();

            // if we are the lower right hand corner of the table,
            // help out the user by moving focus
            if( row+1 == num_rows && col+1 == num_cols )
            {
                // need to complete editing if necessary:
                table.getCellEditor().stopCellEditing();
                
                // move cursor to the first text box:
                next_comp.requestFocusInWindow();
                // deselect cell of table:
                table.getSelectionModel().clearSelection();
                return;
            }
            
            // if we are not at lower right, then invoke the original action:
            invokeOriginalAction( e );

            //  Keep invoking the original action until we find an editable cell
            while (! table.isCellEditable(row, col) )
            {
                invokeOriginalAction( e );


                //  We didn't move anywhere, reset cell selection and get out.

                if (row == table.getSelectedRow()
                &&  col == table.getSelectedColumn())
                {
                        table.changeSelection(originalRow, originalColumn, false, false);
                        break;
                }

                row = table.getSelectedRow();
                col = table.getSelectedColumn();

                //  Back to where we started, get out.

                if (row == originalRow
                &&  col == originalColumn)
                {
                        break;
                }
            }
        }
}