param(
  [Parameter(ValueFromRemainingArguments=$true)]
  [string[]]$Args
)

# Run without TTY (-t) flag to avoid "input device is not a TTY" error
docker run --rm -v ${PWD}:/workspace sage-cursor sage @Args
