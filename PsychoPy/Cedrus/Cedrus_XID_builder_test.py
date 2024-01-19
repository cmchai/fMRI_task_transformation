from psychopy import core, visual, logging, event, parallel, gui
from psychopy.hardware import keyboard
import time
import sys

# Create a keyboard (e.g. to check for escape)
keyB = keyboard.Keyboard()

# Set the log level
logging.console.setLevel(logging.INFO)  # this outputs to the screen, not a file

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
sys.stdout.flush()

# Clear buffer
while len(cedrusBox.response_queue):
    cedrusBox.clear_response_queue()
    cedrusBox.poll_for_response() #often there are more resps waiting!

print('Device ready, press a key')
sys.stdout.flush()

cntr = 0;

while True:
    # Check for responses
    cedrusBox.poll_for_response()   
    
    # Process the responses
    while len(cedrusBox.response_queue):
        
        # Get the value of the next response in the queue
        evt = cedrusBox.get_next_response()
        
        # Ignore a response not in the list
        if evt['key'] not in [0, 1, 2, 3, 4]:
            continue  # we don't care about this key
            
        # Only process pressed keys, not released keys
        if evt['pressed']:
            print('Event ' + str(cntr) + ': ' + str(evt['key']))
            cntr += 1
            
            # Empty the response queue and stop processing
            cedrusBox.poll_for_response()
            cedrusBox.clear_response_queue()  # don't process again

    # Quit the loop if escape has been pushed
    if keyB.getKeys(keyList=["escape"]):
        core.quit()