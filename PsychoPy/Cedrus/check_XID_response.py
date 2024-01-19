from psychopy import core, visual, logging, event, parallel, gui
from psychopy.hardware import keyboard

def check_XID_response(cedrusBox, maxWait, keyList, responseClock):
    responseClock.reset()
    respKey = 4 # if no response from participant all, the only response will be the trigger 
    # clear buffer
    while len(cedrusBox.response_queue):
        cedrusBox.clear_response_queue()
        cedrusBox.poll_for_response() #often there are more resps waiting!

    while responseClock.getTime() <= maxWait: # maxWait
        # Check for responses
        cedrusBox.poll_for_response()

        # Process the responses
        if cedrusBox.has_response():
            respTime = responseClock.getTime();
            
            # Get the value of the next response in the queue
            evt = cedrusBox.get_next_response()

            # Ignore a response not in the list
            if evt['key'] not in keyList:
                continue  # we don't care about this key
            
            # Only process pressed keys, not released keys
            if evt['pressed']:
                respKey = evt['key']
                # Empty the response queue and stop processing
                cedrusBox.poll_for_response()
                cedrusBox.clear_response_queue()  # don't process again
                break

    return respKey, respTime
    
