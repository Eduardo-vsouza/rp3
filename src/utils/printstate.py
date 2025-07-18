def print_state(message, level=1, color='default', marker='-'):
    markers = {1: marker*1, 
                2: marker*2,
                3: marker*3,
                4: marker*4,
        }
    colors = {
        'green': '\033[92m',
        'red': '\033[91m',
        'yellow': '\033[93m',
        'blue': '\033[94m',
    }
    if color == 'default':
        message = f'{markers[level]} {message}'
    else:
        message = f'{colors[color]}{markers[level]} {message}\033[0m'
    print(message)