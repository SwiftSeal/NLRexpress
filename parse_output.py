import sys

def process_input(input_data):
    # Implement your logic to process the input data here
    # For demonstration purposes, let's just print the input
    print("Input received:", input_data)

def main():
    try:
        # Read input from stdin
        input_data = sys.stdin.read().strip()

        # Process the input data
        process_input(input_data)

    except KeyboardInterrupt:
        # Handling KeyboardInterrupt (Ctrl+C) gracefully
        print("Script terminated by user.")
        sys.exit(0)

if __name__ == "__main__":
    main()
