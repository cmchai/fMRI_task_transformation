from psychopy import core, visual, logging, event, parallel, gui
from psychopy.hardware import keyboard
import time

# Create a keyboard (e.g. to check for escape)
keyB = keyboard.Keyboard()

# Create a Cedrus device
# code copied from Builder file
try: 
   import pyxid2 as pyxid
except ImportError:
   import pyxid
cedrusBox = None

for n in range(10):  # doesn't always work first time!
    try:
        devices = pyxid.get_xid_devices()
        core.wait(0.1)
        cedrusBox = devices[0]
        break  # found the device so can break the loop
    except Exception:
        pass
if not cedrusBox:
    logging.error('Could not find a Cedrus device.')
    core.quit()

print('Cedrus device found: ' + cedrusBox.device_name)

# Clear buffer
while len(cedrusBox.response_queue):
    cedrusBox.clear_response_queue()
    cedrusBox.poll_for_response() #often there are more resps waiting!

cntr = 0;

cedrusBox.reset_rt_timer()
cedrusBox.clear_response_queue()
while True:
    # Check for responses
    cedrusBox.poll_for_response()

    # Process the responses
    if cedrusBox.has_response():
        # Get the value of the next response in the queu
        evt = cedrusBox.get_next_response()

        # Do some housekeeping
        if evt['pressed']:
            evtType = 'Pressed'
        else:
            evtType = 'Released'

        print('Event ' + str(cntr) + ' | ' + 
              'Key: ' + str(evt['key']) + ' | ' + 
              'Type: ' + evtType + ' | ' + 
              'Time: ' + str(evt['time']))
        cntr += 1

    # Quit the loop if escape has been pushed
    if keyB.getKeys(keyList=["escape"]):
        core.quit()
