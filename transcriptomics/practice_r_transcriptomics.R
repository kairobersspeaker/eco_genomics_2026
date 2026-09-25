## Playing in R to understand Basic R functions and Objects####

x <- 5
students <- data.frame(
  name = c("A","B","C"),
  height = c(62,68,72)
)
head(students)
class(students)
str(students)

students$name[2]

mean(students$height)


