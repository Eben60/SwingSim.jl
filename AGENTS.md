# ShareAdd.jl - AI Agent Technical Reference

## Agents Behavior

- **Always clarify first** if a user's request is unclear, before starting the actual action.
- **Do not** interpret a question or a review request as an implicit request for action. Example of proper dialogue:
    - *Human*: Is XY a good idea?
    - *Agent*: Yes, XY is good because of A, B, and C. Should I implement it for you?
    - *Human*: Yes, please
- **Do not** update this file unless explicitely requested.
- Use Kaimon MCP where expedient. If it is not available, while you may need it, pause and let me know
- If I tell you "aopp", it means: "Ask (if you have any questions) Otherwise Please Proceed"
