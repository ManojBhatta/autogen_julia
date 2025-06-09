from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
from autogen.retrieve_utils import LocalCodebaseRetreiver 
from autogen_julia.julia_exec.julia_exec import create_julia_agents
from autogen import GroupChat, GroupChatManager

retriever = LocalCodebaseRetriever(code_dir="/home/manoj/ANK/autogen/autogen_julia/TemperFEM_rag_inp")
rag_agent = RetrieveUserProxyAgent(
    name="rag_agent",
    system_message="You are a helpful assistant that understands code.",
    code_retriever=retriever
)

user_proxy, _ = create_julia_agents()

groupchat = GroupChat(agents=[rag_agent, user_proxy], messages=[], max_round=10)

manager = GroupChatManager(groupchat=groupchat)

while True:
    user_input = input("\nWhat Julia code would you like to create and run? ")
    
    if user_input.lower() in ['quit', 'exit', 'q']:
        break
    
    try:
        # Use the coder to generate code, then execute it
        user_proxy.initiate_chat(
            rag_agent,
            message=user_input,
            max_turns=10
        )
        
    except KeyboardInterrupt:
        print("\nDemo interrupted by user.")
        break
    except Exception as e:
        print(f"Error: {e}")