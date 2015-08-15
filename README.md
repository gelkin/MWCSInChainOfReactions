Course work. Find mwcs in chain of reactions
=======================
# compile
```python
mvn compile
```
# run
```python
mvn exec:java -Dexec.args="nodes.txt edges.txt signals.txt initial_solution.txt"
```
> -- last argument is optional

# result
After running results will be at:
- src/main/resources/result-nodes.txt
- src/main/resources/result-edges.txt
