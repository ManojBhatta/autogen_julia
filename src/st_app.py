import streamlit as st
from julia_exec import create_julia_agents
import os
import time

st.title("Simulation Agents")

user, julia_coder = create_julia_agents()

if "model" not in st.session_state:
    st.session_state["model"] = julia_coder
if "user" not in st.session_state:
    st.session_state["user"] = user
if "messages" not in st.session_state:
    st.session_state["messages"] = []
if "last_files" not in st.session_state:
    st.session_state["last_files"] = set(os.listdir("results"))

def get_new_files():
    """Detect new image/video files in results folder."""
    exts = (".png", ".jpg", ".jpeg", ".gif", ".mp4", ".webm")
    all_files = set(f for f in os.listdir("results") if f.lower().endswith(exts))
    new_files = all_files - st.session_state["last_files"]
    st.session_state["last_files"] = all_files
    return sorted(new_files)

def is_execution_successful(message):
    """Check if execution was successful by looking for exit_code: 0"""
    return "exit_code=0" in message or "exitcode: 0" in message

# Display previous messages and any associated images/animations
for message in st.session_state["messages"]:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])
        # If message has associated files, display them
        if "files" in message:
            for file in message["files"]:
                file_path = os.path.join("results", file)
                if file.lower().endswith((".png", ".jpg", ".jpeg", ".gif")):
                    st.image(file_path)
                elif file.lower().endswith((".mp4", ".webm")):
                    st.video(file_path)

if prompt := st.chat_input("Enter your query"):
    st.session_state["messages"].append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    client = st.session_state["model"]
    user = st.session_state["user"]
    
    # Build conversation from all previous messages for context
    conversation = []
    for msg in st.session_state["messages"]:
        conversation.append({"role": msg["role"], "content": msg["content"]})
    
    max_iterations = 5  # Prevent infinite loops
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        
        # Get code response from the Julia agent
        response = client.generate_reply(conversation, sender=user)
        
        with st.chat_message("assistant"):
            st.markdown(response)
        
        st.session_state["messages"].append({
            "role": "assistant", 
            "content": response,
            "files": []
        })
        
        # Check if agent wants to terminate
        if "TERMINATE" in response:
            break
            
        # Add response to conversation
        conversation.append({"role": "assistant", "content": response})
        
        # Execute the code
        exec_response = user.generate_reply(conversation, sender=client)
        
        with st.chat_message("user"):
            st.markdown(exec_response)
            
            # Wait for files to be written
            time.sleep(1)
            
            # Detect new files generated in results
            new_files = get_new_files()
            for file in new_files:
                file_path = os.path.join("results", file)
                if file.lower().endswith((".png", ".jpg", ".jpeg", ".gif")):
                    st.image(file_path)
                elif file.lower().endswith((".mp4", ".webm")):
                    st.video(file_path)
        
        # Save execution result and files
        st.session_state["messages"].append({
            "role": "user",
            "content": exec_response,
            "files": list(new_files) if new_files else []
        })
        
        # Add execution result to conversation
        conversation.append({"role": "user", "content": exec_response})
        
        # Check if execution was successful
        if is_execution_successful(exec_response):
            break