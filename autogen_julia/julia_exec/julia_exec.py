import subprocess
import os
import tempfile
import re
from typing import List, Optional, Union
from pathlib import Path
import autogen
from autogen.coding import CodeBlock, CodeExecutor, CodeResult, CodeExtractor
from autogen import config_list_from_json

script_dir = Path(__file__).resolve().parent
llm_config = config_list_from_json(f"{script_dir}/config.json")

class JuliaCodeExtractor:
    """Code extractor for Julia code blocks"""
    
    def extract_code_blocks(self, message: str | list | None) -> List[CodeBlock]:
        """Extract Julia code blocks from a message"""
        if message is None:
            return []
        
        # Handle list input (for multimodal messages)
        if isinstance(message, list):
            # Extract text from message parts
            text_content = ""
            for part in message:
                if hasattr(part, 'text'):
                    text_content += part.text + "\n"
                elif isinstance(part, dict) and 'text' in part:
                    text_content += part['text'] + "\n"
            message = text_content
        
        if not isinstance(message, str):
            return []
        
        # Pattern to match code blocks with language specification
        pattern = r'```(?:julia|jl)\n(.*?)```'
        matches = re.findall(pattern, message, re.DOTALL)
        
        code_blocks = []
        for match in matches:
            code_blocks.append(CodeBlock(code=match.strip(), language="julia"))
        
        # Also check for generic code blocks if no julia-specific blocks found
        if not code_blocks:
            generic_pattern = r'```\n(.*?)```'
            generic_matches = re.findall(generic_pattern, message, re.DOTALL)
            for match in generic_matches:
                # Simple heuristic: if it contains Julia-like syntax, assume it's Julia
                if any(keyword in match for keyword in ['println', 'function', 'end', 'using']):
                    code_blocks.append(CodeBlock(code=match.strip(), language="julia"))
        
        return code_blocks


class JuliaCodeExecutor(CodeExecutor):
    """Custom Julia Code Executor for AutoGen"""
    
    def __init__(self, 
                 timeout: int = 60,
                 work_dir: Optional[Union[Path, str]] = "results",
                 julia_executable: str = "julia"):
        """
        Initialize Julia Code Executor
        Args:
            timeout: Maximum execution time in seconds
            work_dir: Working directory for code execution
            julia_executable: Path to Julia executable
        """
        self.timeout = timeout
        self.julia_executable = julia_executable
        self._code_extractor = JuliaCodeExtractor()
        
        if work_dir is None:
            self.work_dir = Path(tempfile.mkdtemp())
        else:
            self.work_dir = Path(work_dir)
            
        self.work_dir.mkdir(exist_ok=True, parents=True)
        
        # Check if Julia is available
        self._check_julia_installation()
    
    @property
    def code_extractor(self) -> CodeExtractor:
        """Return the code extractor for this executor"""
        return self._code_extractor
    
    def _check_julia_installation(self):
        """Check if Julia is properly installed"""
        try:
            result = subprocess.run([self.julia_executable, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                raise RuntimeError(f"Julia not found or not working: {result.stderr}")
            print(f"Julia version detected: {result.stdout.strip()}")
        except FileNotFoundError:
            raise RuntimeError("Julia executable not found. Please install Julia and ensure it's in your PATH.")
        except subprocess.TimeoutExpired:
            raise RuntimeError("Julia installation check timed out.")
    
    def execute_code_blocks(self, code_blocks: List[CodeBlock]) -> CodeResult:
        """Execute Julia code blocks"""
        
        if not code_blocks:
            return CodeResult(exit_code=0, output="No code to execute")
        
        outputs = []
        exit_code = 0
        
        for i, code_block in enumerate(code_blocks):
            if code_block.language.lower() not in ["julia", "jl"]:
                outputs.append(f"Skipping non-Julia code block (language: {code_block.language})")
                continue
            
            try:
                result = self._execute_julia_code(code_block.code, f"block_{i}")
                outputs.append(f"=== Code Block {i+1} ===")
                outputs.append(result.output)
                
                if result.exit_code != 0:
                    exit_code = result.exit_code
                    
            except Exception as e:
                outputs.append(f"Error executing code block {i+1}: {str(e)}")
                exit_code = 1
        
        return CodeResult(
            exit_code=exit_code,
            output="\n".join(outputs)
        )
    
    def _execute_julia_code(self, code: str, filename_prefix: str = "temp") -> CodeResult:
        """Execute a single Julia code block"""
        
        # Create temporary Julia file
        # julia_file =  f"{filename_prefix}.jl"
        julia_file = os.path.join(self.work_dir, f"{filename_prefix}.jl")
        print(f'\n \n the path to julia file is {julia_file} \n \n ')
        
        try:
            # Write code to file
            with open(julia_file, 'w', encoding='utf-8') as f:
                f.write(code)
            
            # Execute Julia code
            cmd = [self.julia_executable, str(f"{filename_prefix}.jl")]
            
            result = subprocess.run(
                cmd,
                cwd=self.work_dir,
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            # Combine stdout and stderr
            output = ""
            if result.stdout:
                output += f"Output:\n{result.stdout}"
            if result.stderr:
                output += f"\nErrors/Warnings:\n{result.stderr}"
            
            return CodeResult(
                exit_code=result.returncode,
                output=output if output else "Code executed successfully (no output)"
            )
            
        except subprocess.TimeoutExpired:
            return CodeResult(
                exit_code=124,
                output=f"Julia code execution timed out after {self.timeout} seconds"
            )
        except Exception as e:
            return CodeResult(
                exit_code=1,
                output=f"Execution error: {str(e)}"
            )
        # finally:
        #     # Clean up temporary file
        #     if julia_file.exists():
        #         julia_file.unlink()
    
    def restart(self) -> None:
        """Restart the code executor (cleanup temporary files)"""
        # Clean up any remaining temporary files
        for file in self.work_dir.glob("*.jl"):
            try:
                file.unlink()
            except:
                pass


def create_julia_agents():
    """Create AutoGen agents for Julia development"""
    
    # Create Julia Code Executor
    julia_executor = JuliaCodeExecutor(timeout=600)
    
    user_proxy = autogen.UserProxyAgent(
        name="user_proxy",
        human_input_mode="TERMINATE",
        max_consecutive_auto_reply=10,
        is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config={
            "executor": julia_executor,
        },
        # system_message=f"you are a user proxy agent, you act as a virtual version of the user and check the execution result of the code  provided by julia coder and executed by the code executor agent, if the code is ran successfully, just send the message 'TERMINATE' and nothing else and if any error is there suggest what might be causing the error and possible modificatoins  "
    )

#     julia_executor_agent = autogen.UserProxyAgent(
#     name="julia_executor_agent",
#     human_input_mode="NEVER",  
#     # max_consecutive_auto_reply=3,
#     is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
#     code_execution_config={
#         "executor": julia_executor,
#     },
#     system_message=f"You automatically execute any Julia code that is provided to you and report the results. and if any plots are to be generated or any file is to be saved , save the file with appropraite name to {julia_executor.work_dir}"
# )
# This is redundant as the user proxy can execute the code
   
    # Julia Code Generator Agent
    julia_coder = autogen.AssistantAgent(
        name="julia_coder",
        llm_config=llm_config,
        system_message="""You are a Julia programming expert. Your role is to:
        1. Generate high-quality Julia code based on user requirements
        2. Explain Julia concepts and best practices
        3. Debug and optimize Julia code
        4. Suggest appropriate Julia packages for specific tasks
        
        When writing Julia code:
        - Use clear, readable syntax
        - Include appropriate comments
        - Follow Julia naming conventions
        - Use type annotations where helpful
        - Consider performance implications
        
        Always wrap your Julia code in ```julia blocks.
        the following is only when you receive the execution result from the user_proxy agent, you should not execute the code yourself, you just generate the code and send it to the user_proxy agent for execution,
        and if the code is executed successfully, you should send the message 'TERMINATE' and nothing else,
        and if there is any error in the code, you should suggest what might be causing the error and possible modifications to fix it,
        """
    )
    
    return user_proxy, julia_coder


# Example usage function
def main():
    """Main function to demonstrate the Julia AutoGen system"""
    
    # Create agents
    user_proxy, julia_coder = create_julia_agents()
    
    # Start the conversation
    print("\n=== Julia AutoGen System Started ===")
    print("You can now interact with the Julia development team!")
    print("Available agents:")
    print("- julia_coder: Generates Julia code")
    print("- julia_executor: Executes and analyzes Julia code")
    print("\nType 'TERMINATE' to end the conversation.")
    
    # Example group chat
    groupchat = autogen.GroupChat(
        agents=[user_proxy, julia_coder],
        messages=[],
        max_round=50
    )
    
    manager = autogen.GroupChatManager(groupchat=groupchat, llm_config=llm_config)
    
    # Start the conversation
    user_proxy.initiate_chat(
        manager,
        message="Hello! I need help with Julia programming. What can you help me with?"
    )


if __name__ == "__main__":
    main()