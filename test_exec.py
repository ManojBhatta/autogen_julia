#!/usr/bin/env python3
"""
Julia AutoGen Workflow Demo
This script demonstrates the proper workflow for generating and executing Julia code
"""

from autogen_julia.julia_exec.julia_exec import create_julia_agents
from pathlib import Path
from autogen import config_list_from_json

script_dir = Path(__file__).resolve().parent
llm_config = config_list_from_json(f"{script_dir}/config.json")

def demo_interactive_session():
    """Demo interactive session with user input"""
    print("\n💬 Demo: Interactive Session")
    print("=" * 50)
    
    # Create agents
    user_proxy, julia_coder = create_julia_agents()
    
    print("\nStarting interactive session...")
    print("You can ask for Julia code, and it will be generated and executed automatically.")
    print("Type 'quit' to exit this demo.")
    
    while True:
        user_input = input("\nWhat Julia code would you like to create and run? ")
        
        if user_input.lower() in ['quit', 'exit', 'q']:
            break
        
        try:
            # Use the coder to generate code, then execute it
            user_proxy.initiate_chat(
                julia_coder,
                message=user_input,
                max_turns=10
            )
            
        except KeyboardInterrupt:
            print("\nDemo interrupted by user.")
            break
        except Exception as e:
            print(f"Error: {e}")
    
    print("Interactive demo ended.")

def main():
    """Main demo function"""
    print("🎯 Julia AutoGen System Demo")
    print("=" * 60)
    
    try:
        demo_interactive_session()
        
    except Exception as e:
        print(f"❌ Demo failed: {e}")

if __name__ == "__main__":
    main()